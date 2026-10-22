# Float / oxygen_saturation.py

Questo README descrive il ruolo e il funzionamento di `oxygen_saturation.py` nella cartella `bit.sea/src/bitsea/Float`.

## Cosa fa `oxygen_saturation.py`

`oxygen_saturation.py` contiene due funzioni principali:

1. `oxy_sat(p)`
   - calcola la saturazione dell'ossigeno in acqua (`O2o` in mmol/m3)
   - usa la formula empirica di Garcia e Gordon (1992)
   - calcola la temperatura media e la salinità media nella colonna superiore fino a 5 m di profondità
   - se non trova sufficiente dati di `TEMP` o `PSAL` nella frazione superficiale, restituisce `nan`

2. `oxy_check(Pres, Prof, p)`
   - confronta il valore di ossigeno superficiale misurato (`Osurf`) con il valore di saturazione teorico (`Osat`)
   - definisce `Osurf` come la media dei valori `DOXY` a pressione <= 5 m
   - considera il profilo valido solo se la differenza assoluta è minore di 10 mmol/m3
   - restituisce `True` se il profilo passa il controllo, `False` in tutti gli altri casi

Questa verifica è usata come test di qualità sui profili `DOXY` per evitare di salvare profili con valori di ossigeno superficiale fisicamente improbabili.

## Workflow di utilizzo in `superfloat_oxygen.py`

Il controllo di saturazione viene usato durante il caricamento dei profili CORIOLIS e prima di scrivere un record in SUPERFLOAT.

Diagramma semplificato del workflow:

```
START
  |
  v
Carica lista profili CORIOLIS
  |
  v
Per ogni profilo:
  - leggi DOXY
  - leggi TEMP e PSAL
  - chiama oxy_check()
      |\
      | `--> False -> rigetta profilo
      |
      `--> True -> conserva profilo
  |
  v
Se il profilo è valido:
  - lo inserisci nel dataset di lavoro
  - può essere successivamente usato per il drift check
  - può essere scritto in SUPERFLOAT se passa gli altri controlli
```

## Dettaglio del controllo

- `oxy_sat(p)` calcola la saturazione teorica
  - input: profilo `p` con variabili `TEMP` e `PSAL`
  - output: saturazione ossigeno in mmol/m3
- `oxy_check(Pres, Prof, p)` calcola il valore superficiale misurato
  - input:
    - `Pres`: vettore delle pressioni del profilo DOXY
    - `Prof`: vettore dei valori DOXY
    - `p`: oggetto profilo per leggere `TEMP` e `PSAL`
  - output: `True` / `False`

Condizioni di rifiuto:

- nessun valore `TEMP` o `PSAL` entro 5 m
- nessun valore `DOXY` entro 5 m
- differenza tra `Osurf` e saturazione > 10 mmol/m3

## Prima vs Ora

### Prima
- il controllo sulla saturazione poteva essere gestito direttamente all'interno del codice principale di generazione SUPERFLOAT
- la logica era meno modulare e difficile da riusare per altri script
- la verifica era parte di un flusso monolitico, più difficile da isolare e testare

### Ora
- `oxygen_saturation.py` è una utility separata per il controllo di qualità di `DOXY`
- il workflow di `superfloat_oxygen.py` la utilizza come filtro standard su profili aggiustati
- questo rende il comportamento più chiaro, riutilizzabile e manutentibile
- le eccezioni vengono gestite con un unico punto di controllo e il risultato `True/False` determina l'accettazione del profilo

## Perché è importante

Questo modulo è essenziale per evitare che profili con ossigeno superficiale anomalo vengano inclusi nel dataset SUPERFLOAT. Il test di saturazione è una barriera fisica: un valore DOXY molto diverso dalla saturazione teorica indica probabilmente un errore di misura o un profilo non affidabile.

## Nota

Il file `oxygen_saturation.py` non gestisce direttamente la scrittura dei file; si limita a valutare la qualità del profilo. La scrittura viene fatta da `superfloat_oxygen.py` in base al risultato del controllo.
