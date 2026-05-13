#!/usr/bin/python3
import csv
import io
import re
import zipfile
import unicodedata
import warnings
from pathlib import Path

import pandas as pd
import requests
import geopandas as gpd


# =========================================================
# Konfigurace
# =========================================================
OUTPUT_FILE = "mesta_cr.txt"
WORKDIR = Path("mesta_cr_cache")
WORKDIR.mkdir(exist_ok=True)

# ČÚZK / RÚIAN
UI_OBEC_ZIP_URL = "https://services.cuzk.cz/sestavy/cis/UI_OBEC.zip"
CS_STATUS_OBCE_ZIP_URL = "https://services.cuzk.cz/sestavy/cis/CS_STATUS_OBCE.zip"
RUIAN_STATE_SHP_ZIP_URL = "https://services.cuzk.gov.cz/shp/stat/epsg-5514/1.zip"

# ČSÚ MOS open data (oficiální otevřená data za obce)
MOS_DATA_2024_URL = "https://opendata.csu.gov.cz/soubory/od/od_mos01/mos_data_2024.csv"
MOS_UKAZ_URL = "https://opendata.csu.gov.cz/soubory/od/od_mos01/mos_ukaz.csv"

HEADERS = {
    "User-Agent": "mesta-cr-export/2.1"
}


# =========================================================
# Pomocné funkce
# =========================================================
def download_bytes(url: str) -> bytes:
    r = requests.get(url, headers=HEADERS, timeout=180)
    r.raise_for_status()
    return r.content


def download_file(url: str, path: Path) -> None:
    """
    Stáhne soubor jen pokud ještě neexistuje.
    U ZIPu navíc ověří, že začíná PK magic bytes.
    """
    if path.exists() and path.stat().st_size > 0:
        return

    data = download_bytes(url)

    if url.lower().endswith(".zip"):
        if len(data) < 4 or data[:2] != b"PK":
            preview = data[:200].decode("utf-8", errors="replace")
            raise RuntimeError(
                f"Stažený soubor z {url} nevypadá jako ZIP archiv. "
                f"Začátek obsahu: {preview!r}"
            )

    path.write_bytes(data)


def normalize(text) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(text).strip().lower())


def normalize_name(text: str) -> str:
    s = str(text).strip().lower()
    s = unicodedata.normalize("NFKD", s)
    s = "".join(ch for ch in s if not unicodedata.combining(ch))

    prefixes = [
        "hlavni mesto ",
        "hl m ",
        "hl m",
        "mesto ",
        "mestys ",
        "obec ",
        "statutarni mesto ",
    ]
    for p in prefixes:
        if s.startswith(p):
            s = s[len(p):]

    s = s.replace("–", "-").replace("—", "-")
    s = re.sub(r"\(.*?\)", "", s)
    s = re.sub(r"[^a-z0-9]+", "", s)
    return s


def detect_delimiter(text: str) -> str:
    """
    Pokusí se odhadnout oddělovač.
    Když sniff selže, vrátí čárku jako rozumný default.
    """
    sample = "\n".join(text.splitlines()[:20])
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=";,|\t,")
        return dialect.delimiter
    except Exception:
        return ","


def decode_bytes(raw: bytes) -> str:
    for enc in ("utf-8-sig", "cp1250", "utf-8", "latin1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            pass
    raise RuntimeError("Nepodařilo se dekódovat textový soubor.")


def read_csv_url_to_df(url: str) -> pd.DataFrame:
    """
    Robustní načtení CSV z URL.

    Postup:
    1) stáhne bytes a dekóduje text
    2) zkusí více oddělovačů přes pandas
    3) při chybách přeskočí poškozené řádky
    4) když pandas selže, použije fallback přes csv.reader(strict=False)
    """
    raw = download_bytes(url)
    text = decode_bytes(raw)

    candidate_delimiters = []
    sniffed = detect_delimiter(text)
    if sniffed:
        candidate_delimiters.append(sniffed)

    for d in [",", ";", "\t", "|"]:
        if d not in candidate_delimiters:
            candidate_delimiters.append(d)

    last_error = None

    # 1) pokus přes pandas
    for delimiter in candidate_delimiters:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                df = pd.read_csv(
                    io.StringIO(text),
                    sep=delimiter,
                    dtype=str,
                    engine="python",
                    on_bad_lines="skip",
                    quotechar='"',
                    doublequote=True,
                )
            if df is not None and len(df.columns) >= 2 and len(df) > 0:
                return df
        except Exception as e:
            last_error = e

    # 2) fallback přes csv modul
    for delimiter in candidate_delimiters:
        try:
            rows = []
            reader = csv.reader(
                io.StringIO(text),
                delimiter=delimiter,
                quotechar='"',
                doublequote=True,
                strict=False,
            )

            header = None
            max_cols = 0

            for row in reader:
                if not row:
                    continue
                if header is None:
                    header = row
                    max_cols = len(header)
                    continue

                if len(row) < max_cols:
                    row = row + [""] * (max_cols - len(row))
                elif len(row) > max_cols:
                    row = row[:max_cols]

                rows.append(row)

            if header and rows:
                df = pd.DataFrame(rows, columns=header)
                if len(df.columns) >= 2 and len(df) > 0:
                    return df

        except Exception as e:
            last_error = e

    raise RuntimeError(f"Nepodařilo se naparsovat CSV z {url}. Poslední chyba: {last_error}")


def read_cuzk_csv_from_zip(url: str) -> list[dict]:
    raw_zip = download_bytes(url)

    with zipfile.ZipFile(io.BytesIO(raw_zip)) as zf:
        csv_names = [n for n in zf.namelist() if n.lower().endswith(".csv")]
        if not csv_names:
            raise RuntimeError(f"V ZIPu z {url} nebyl nalezen CSV soubor.")
        raw = zf.read(csv_names[0])

    text = None
    for enc in ("utf-8-sig", "cp1250", "utf-8"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            pass

    if text is None:
        raise RuntimeError(f"Nepodařilo se dekódovat CSV z {url}.")

    reader = csv.DictReader(io.StringIO(text), delimiter=";")
    rows = list(reader)

    if not rows:
        raise RuntimeError(f"CSV z {url} neobsahuje žádné řádky.")

    return rows


def find_col(columns, *needles):
    norm = {c: normalize(c) for c in columns}
    for col, n in norm.items():
        if all(nd in n for nd in needles):
            return col
    return None


# =========================================================
# ČSÚ MOS populace (2024)
# =========================================================
def load_mos_population() -> pd.DataFrame:
    """
    Načte populaci z ČSÚ MOS open data 2024.

    Očekávaná struktura:
      rok, kodukaz, koduzemi, hodnota

    Postup:
    1) načti číselník ukazatelů (mos_ukaz.csv)
    2) najdi nejpravděpodobnější ukazatel pro počet obyvatel
    3) načti mos_data_2024.csv
    4) vyfiltruj daný kodukaz a vrať KOD + POPULATION
    """
    ukaz_df = read_csv_url_to_df(MOS_UKAZ_URL)
    ukaz_df.columns = [str(c).strip() for c in ukaz_df.columns]

    print("Sloupce mos_ukaz.csv:", list(ukaz_df.columns))

    code_col = (
        find_col(ukaz_df.columns, "kod", "ukaz")
        or find_col(ukaz_df.columns, "kod")
    )
    text_cols = [c for c in ukaz_df.columns if c != code_col]

    if code_col is None:
        raise RuntimeError(
            f"V mos_ukaz.csv se nepodařilo najít sloupec s kódem ukazatele. Sloupce: {list(ukaz_df.columns)}"
        )

    tmp = ukaz_df.copy()
    tmp["_SEARCH_TEXT"] = tmp[text_cols].fillna("").astype(str).agg(" | ".join, axis=1)

    candidates = tmp[
        tmp["_SEARCH_TEXT"].str.contains(r"počet obyvatel|pocet obyvatel", case=False, regex=True, na=False)
    ].copy()

    if candidates.empty:
        raise RuntimeError("V mos_ukaz.csv se nepodařilo najít kandidátní ukazatel pro počet obyvatel.")

    bad_patterns = r"narozen|zemřel|zemrel|sňatk|snatk|rozvod|potrat|stěhov|stehov|přistěh|vystěh|průměrný věk|prumerny vek"
    candidates["_BAD"] = candidates["_SEARCH_TEXT"].str.contains(bad_patterns, case=False, regex=True, na=False)

    candidates["_GOOD_SCORE"] = 0
    candidates.loc[
        candidates["_SEARCH_TEXT"].str.contains(r"stav obyvatel", case=False, regex=True, na=False),
        "_GOOD_SCORE"
    ] += 3
    candidates.loc[
        candidates["_SEARCH_TEXT"].str.contains(r"počet obyvatel|pocet obyvatel", case=False, regex=True, na=False),
        "_GOOD_SCORE"
    ] += 2
    candidates.loc[
        candidates["_SEARCH_TEXT"].str.contains(r"celkem", case=False, regex=True, na=False),
        "_GOOD_SCORE"
    ] += 1

    candidates = candidates.sort_values(by=["_BAD", "_GOOD_SCORE"], ascending=[True, False])

    chosen_row = candidates.iloc[0]
    indicator_code = str(chosen_row[code_col]).strip()

    print(f"Použitý ukazatel populace z MOS: {indicator_code}")
    print(f"Text ukazatele: {chosen_row['_SEARCH_TEXT']}")

    data_df = read_csv_url_to_df(MOS_DATA_2024_URL)
    data_df.columns = [str(c).strip() for c in data_df.columns]

    print("Sloupce mos_data_2024.csv:", list(data_df.columns))

    required = {"rok", "kodukaz", "koduzemi", "hodnota"}
    if not required.issubset(set(data_df.columns)):
        raise RuntimeError(
            f"mos_data_2024.csv nemá očekávané sloupce {required}. Skutečné sloupce: {list(data_df.columns)}"
        )

    pop_df = data_df[data_df["kodukaz"].astype(str).str.strip() == indicator_code].copy()

    if pop_df.empty:
        raise RuntimeError(f"V mos_data_2024.csv nebyla nalezena žádná data pro ukazatel {indicator_code}.")

    pop_df["KOD"] = (
        pop_df["koduzemi"]
        .astype(str)
        .str.extract(r"(\d+)", expand=False)
        .fillna("")
        .str.zfill(6)
    )

    pop_df["POPULATION"] = (
        pop_df["hodnota"]
        .astype(str)
        .str.replace("\xa0", "", regex=False)
        .str.replace(" ", "", regex=False)
        .str.extract(r"(\d+)", expand=False)
    )

    pop_df = pop_df.dropna(subset=["KOD", "POPULATION"])
    pop_df = pop_df[pop_df["POPULATION"].str.match(r"^\d+$", na=False)]
    pop_df["POPULATION"] = pop_df["POPULATION"].astype(int)

    pop_df = (
        pop_df[["KOD", "POPULATION"]]
        .groupby("KOD", as_index=False)["POPULATION"]
        .max()
    )

    print(f"Nalezeno populačních řádků z MOS: {len(pop_df)}")
    return pop_df


# =========================================================
# ČÚZK SHP – souřadnice
# =========================================================
def extract_zip(zip_path: Path, extract_dir: Path) -> None:
    if extract_dir.exists() and any(extract_dir.iterdir()):
        return

    extract_dir.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_dir)


def find_candidate_shp_files(root_dir: Path) -> list[Path]:
    return sorted(root_dir.rglob("*.shp"))


def choose_best_municipality_shp(shp_files: list[Path]) -> Path:
    last_error = None

    prioritized = sorted(
        shp_files,
        key=lambda p: (
            0 if re.search(r"obec|obce|obc", p.name.lower()) else 1,
            len(str(p))
        )
    )

    for shp in prioritized:
        try:
            gdf = gpd.read_file(shp)
            cols_norm = {c: normalize(c) for c in gdf.columns}

            has_code = any(v == "kod" or v.endswith("kod") for v in cols_norm.values())
            has_name = any("nazev" in v for v in cols_norm.values())
            has_poly = gdf.geometry.geom_type.isin(["Polygon", "MultiPolygon"]).any()

            if has_code and has_name and has_poly:
                return shp

        except Exception as e:
            last_error = e

    raise RuntimeError(
        "Nepodařilo se najít vhodný SHP soubor s polygony obcí."
    ) from last_error


def load_municipality_points_from_shp(zip_path: Path) -> pd.DataFrame:
    extract_dir = zip_path.with_suffix("")
    extract_zip(zip_path, extract_dir)

    shp_files = find_candidate_shp_files(extract_dir)
    if not shp_files:
        raise RuntimeError(
            f"Po rozbalení archivu {zip_path} nebyl nalezen žádný .shp soubor."
        )

    shp_path = choose_best_municipality_shp(shp_files)
    print(f"Použitý SHP soubor: {shp_path}")

    gdf = gpd.read_file(shp_path)

    code_col = None
    for c in gdf.columns:
        n = normalize(c)
        if n == "kod" or n.endswith("kod"):
            code_col = c
            break

    if code_col is None:
        raise RuntimeError(
            f"Ve vrstvě obcí nebyl nalezen sloupec s kódem. Sloupce: {list(gdf.columns)}"
        )

    pts = gdf.copy()
    pts["geometry"] = pts.representative_point()

    if pts.crs is None:
        pts = pts.set_crs(epsg=5514)

    pts = pts.to_crs(epsg=4326)

    out = pd.DataFrame({
        "KOD": pts[code_col].astype(str).str.extract(r"(\d+)", expand=False).fillna("").str.zfill(6),
        "longitude": pts.geometry.x,
        "latitude": pts.geometry.y,
    })

    out = out[out["KOD"] != ""].drop_duplicates(subset=["KOD"])
    return out


# =========================================================
# 1) RÚIAN / ČÚZK: obce + statusy
# =========================================================
obce_rows = read_cuzk_csv_from_zip(UI_OBEC_ZIP_URL)
status_rows = read_cuzk_csv_from_zip(CS_STATUS_OBCE_ZIP_URL)

status_map = {}
for row in status_rows:
    kod = str(row.get("KOD", "")).strip()
    nazev = str(row.get("NAZEV", "")).strip()
    if kod:
        status_map[kod] = nazev

mesta = []
for row in obce_rows:
    if str(row.get("PLATI_DO", "")).strip():
        continue

    kod = str(row.get("KOD", "")).strip().zfill(6)
    nazev = str(row.get("NAZEV", "")).strip()
    status_kod = str(row.get("STATUS_KOD", "")).strip()
    status_name = status_map.get(status_kod, "").strip().lower()

    if "město" in status_name:
        mesta.append({
            "KOD": kod,
            "NAZEV": nazev,
            "STATUS": status_name,
            "NAME_KEY": normalize_name(nazev),
        })

mesta_df = pd.DataFrame(mesta).drop_duplicates(subset=["KOD"]).copy()

if mesta_df.empty:
    raise RuntimeError("Filtr měst vrátil 0 řádků.")


# =========================================================
# 2) ČSÚ: populace (MOS 2024)
# =========================================================
pop_df = load_mos_population()

result = mesta_df.merge(
    pop_df.drop_duplicates(subset=["KOD"]),
    on="KOD",
    how="left"
)


# =========================================================
# 3) ČÚZK SHP: souřadnice
# =========================================================
local_shp_zip = WORKDIR / "ruian_state_epsg5514.zip"
download_file(RUIAN_STATE_SHP_ZIP_URL, local_shp_zip)

coords_df = load_municipality_points_from_shp(local_shp_zip)

result = (
    result
    .merge(coords_df, on="KOD", how="left")
    .sort_values("NAZEV")
    .reset_index(drop=True)
)


# =========================================================
# 4) Výstup
# =========================================================
with open(OUTPUT_FILE, "w", encoding="utf-8", newline="") as f:
    f.write("nazev_mesta|pocet_obyvatel|latitude|longitude\n")
    for _, row in result.iterrows():
        pop = "" if pd.isna(row.get("POPULATION")) else int(row["POPULATION"])
        lat = "" if pd.isna(row.get("latitude")) else f"{row['latitude']:.6f}"
        lon = "" if pd.isna(row.get("longitude")) else f"{row['longitude']:.6f}"
        f.write(f"{row['NAZEV']}|{pop}|{lat}|{lon}\n")

print("Hotovo:", OUTPUT_FILE)
print("Počet měst:", len(result))
print("Chybějící populace:", int(result["POPULATION"].isna().sum()))
print("Chybějící souřadnice:", int(result["latitude"].isna().sum()))