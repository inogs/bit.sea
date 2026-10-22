# Readme_superfloat_oxygen

## Scopo

Questo file descrive il workflow di `superfloat_oxygen.py` nella cartella `bit.sea/src/bitsea/Float`.

L’obiettivo è spiegare come i profili DOXY vengono:
- controllati
- classificati
- scritti in SUPERFLOAT
- registrati in `Floats_accepted.csv` o `Floats_rejected.csv`

---

## Workflow generale

1. Carico i profili CORIOLIS
2. Applico i controlli preliminari a ogni profilo:
   - `RT` (rapporto rigetto)
   - `PresNone` (profili con pochi dati di pressione)
   - `SaturationTest` (controllo saturazione O2 in superficie)
3. Solo i profili che superano questi check entrano nel dataset principale
4. Per ogni float (WMO) calcolo il drift storico
5. Per ogni profilo valido:
   - confronto l’O2 a 600 m con la climatologia EMODNET del sottobacino
   - se manca il dato a 600 m (`NoClimValue`) lo accetto comunque senza detrend
   - se l’offset supera la soglia, lo rigetto
   - altrimenti lo scrivo in SUPERFLOAT

---

## Cosa succede in `superfloat_oxygen.py`

### 1. Caricamento e pre-filtri
- Profilo con `status_var('DOXY') == 'R'` → escluso
- Profilo con `read_doxy(p)` che restituisce `None` → escluso
- Profilo che fallisce `oxygen_saturation.oxy_check(...)` → escluso

Questi profili vanno in:
- `DataMode_and_Saturation_rejection_doxy.csv`

### 2. Calcolo del drift sul float
- `compute_float_trend(...)` utilizza la storia del float
- se il trend è significativo e verificato, può diventare `DRIFT_CODE = 1`
- se è anomalo, `DRIFT_CODE = -2`
- se non è significativo, `DRIFT_CODE = 0`

### 3. Climatology check e scrittura
- `clim_check(...)` confronta il valore O2 a 600 m con EMODNET
- soglia:
  - `threshold = max(2 * stdev, 15)`

Caso 1: `NoClimValue`
- non esiste un valore valido a 600 m nel dataset storico
- con la modifica richiesta:
  - il profilo viene comunque scritto in SUPERFLOAT
  - non viene detrendizzato
  - viene scritto con `DRIFT_CODE = -1`
  - non compare più in `Floats_rejected.csv`

Caso 2: `climatology_offset`
- `abs(OFFSET) >= threshold`
- profilo rigettato
- scritto in `Floats_rejected.csv`

Caso 3: profilo valido
- scritto in `Floats_accepted.csv`
- scritto in file SUPERFLOAT

---

## Cosa contiene ogni CSV

### `DataMode_and_Saturation_rejection_doxy.csv`
- `RT`
- `PresNone`
- `SaturationTest`

### `Floats_rejected.csv`
- solo profili rigettati per climatology offset
- dopo la tua modifica, non contiene più `NoClimValue`

### `Floats_accepted.csv`
- profili accettati e scritti in SUPERFLOAT
- include anche i casi `NoClimValue` con `DRIFT_CODE = -1`

---

## Regole chiave per l’output SUPERFLOAT

- profili rigettati nei pre-check non vengono scritti
- profili con `climatology_offset` non vengono scritti
- profili `NoClimValue` vengono scritti senza detrend
- profili già presenti su disco non vengono riscritti a meno che non usi `--force`

---

## Conclusione

Il tuo workflow è ora:

- filtro preliminare: `RT`, `PresNone`, `SaturationTest`
- drift storico sul float
- controllo climatologia a 600 m
- `NoClimValue` → accettato con `DRIFT_CODE = -1`
- `climatology_offset` → rigettato
- profili validi → accettati e scritti in SUPERFLOAT
