#!/usr/bin/python3

import csv
import io
import time
import zipfile
import requests

# -----------------------------
# Nastavení
# -----------------------------
OUTPUT_FILE = "mesta_cr.txt"

UI_OBEC_ZIP_URL = "https://services.cuzk.cz/sestavy/cis/UI_OBEC.zip"
CS_STATUS_OBCE_ZIP_URL = "https://services.cuzk.cz/sestavy/cis/CS_STATUS_OBCE.zip"

HEADERS = {
    "User-Agent": "mesta-cr-export/1.0 (Python script for educational use)"
}


# -----------------------------
# Pomocné funkce
# -----------------------------
def download_zip_csv_rows(url):
    """Stáhne ZIP z URL, najde v něm první CSV a vrátí list dict řádků."""
    r = requests.get(url, headers=HEADERS, timeout=60)
    r.raise_for_status()

    zf = zipfile.ZipFile(io.BytesIO(r.content))
    csv_names = [name for name in zf.namelist() if name.lower().endswith(".csv")]
    if not csv_names:
        raise RuntimeError(f"V ZIPu z {url} nebyl nalezen žádný CSV soubor.")

    with zf.open(csv_names[0]) as fh:
        raw = fh.read()

    # ČÚZK CSV bývá někdy UTF-8 BOM, jindy legacy kódování; zkusíme obě varianty
    text = None
    for enc in ("utf-8-sig", "cp1250", "utf-8"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue

    if text is None:
        raise RuntimeError(f"Nepodařilo se dekódovat CSV z {url}.")

    # ČÚZK CSV obvykle používá středník
    reader = csv.DictReader(io.StringIO(text), delimiter=";")
    rows = list(reader)

    if not rows:
        raise RuntimeError(f"CSV z {url} bylo načteno, ale neobsahuje žádné řádky.")

    return rows


def geocode_city(name):
    """
    Vrátí (lat, lon) přes Nominatim.
    Hledáme obec/město v ČR.
    """
    params = {
        "q": f"{name}, Česko",
        "format": "jsonv2",
        "limit": 1
    }
    r = requests.get(
        "https://nominatim.openstreetmap.org/search",
        params=params,
        headers=HEADERS,
        timeout=60
    )
    r.raise_for_status()
    data = r.json()

    if not data:
        return "", ""

    return data[0].get("lat", ""), data[0].get("lon", "")


# -----------------------------
# Načtení číselníků
# -----------------------------
obce_rows = download_zip_csv_rows(UI_OBEC_ZIP_URL)
status_rows = download_zip_csv_rows(CS_STATUS_OBCE_ZIP_URL)

# mapa STATUS_KOD -> NAZEV statusu
status_map = {}
for row in status_rows:
    kod = (row.get("KOD") or "").strip()
    nazev = (row.get("NAZEV") or "").strip()
    if kod:
        status_map[kod] = nazev

# -----------------------------
# Filtrace měst
# -----------------------------
mesta = []

for row in obce_rows:
    # přeskočit zaniklé obce
    if (row.get("PLATI_DO") or "").strip():
        continue

    kod = (row.get("KOD") or "").strip()
    nazev = (row.get("NAZEV") or "").strip()
    status_kod = (row.get("STATUS_KOD") or "").strip()
    status_nazev = status_map.get(status_kod, "").strip().lower()

    # chceme opravdu jen města
    if status_nazev in {"město", "statutární město"}:
        mesta.append({
            "kod": kod,
            "nazev": nazev,
            "status": status_nazev
        })

# seřazení podle názvu
mesta.sort(key=lambda x: x["nazev"])

# -----------------------------
# Výstup
# -----------------------------
with open(OUTPUT_FILE, "w", encoding="utf-8", newline="") as f:
    f.write("nazev_mesta|pocet_obyvatel|latitude|longitude\n")

    for i, mesto in enumerate(mesta, start=1):
        name = mesto["nazev"]

        # Počet obyvatel zde stále nemáme z ČSÚ; nechávám placeholder
        population = "NA"

        lat, lon = geocode_city(name)

        f.write(f"{name}|{population}|{lat}|{lon}\n")

        # Nominatim má být volán šetrně
        time.sleep(1)

        print(f"[{i}/{len(mesta)}] {name}")

print("Hotovo:", OUTPUT_FILE)
print("Počet měst:", len(mesta))
