#!/bin/bash
set -e
set -u

# === DEFAULTS ===
N_RUNS=""
VBIAS=""
SAMPLING=""
FOLDER=""
CHANNELS=()
DO_TAU=false
DO_DARK=false
DO_LASER=false
host="aimtti-plh120p-00"
port=9221
N_WAVES=2000

# === ARGUMENT PARSING ===
while [[ $# -gt 0 ]]; do
	case $1 in
		--n_runs) N_RUNS="$2"; shift 2 ;;
		--vbias) VBIAS="$2"; shift 2 ;;
		--sampling) SAMPLING="$2"; shift 2 ;;
		--folder) FOLDER="$2"; shift 2 ;;
		--channels)
			shift
			echo "📥 Parsing channels..."
			while [[ $# -gt 0 && $1 != --* ]]; do
				echo "  ➕ Aggiungo canale: $1"
				CHANNELS+=("$1")
				shift
			done
			;;
		--tau) DO_TAU=true; shift ;;
		--dark) DO_DARK=true; shift ;;
		--laser) DO_LASER=true; shift ;;
		*)
			echo "❌ Unknown option: $1"
			echo "Usage: $0 --n_runs <N> --vbias <V> --sampling <rate> --folder <dir> --channels <list> [--tau|--dark|--laser]"
			exit 1 ;;
	esac
done

# === CHECK ARGOMENTI OBBLIGATORI ===
missing=()
[[ -z "$N_RUNS" ]] && missing+=("--n_runs")
[[ -z "$VBIAS" ]] && missing+=("--vbias")
[[ -z "$SAMPLING" ]] && missing+=("--sampling")
[[ -z "$FOLDER" ]] && missing+=("--folder")
[[ ${#CHANNELS[@]} -eq 0 ]] && missing+=("--channels")

if (( ${#missing[@]} )); then
	echo "❌ Missing required argument(s): ${missing[*]}"
	echo "Usage: $0 --n_runs <N> --vbias <V> --sampling <rate> --folder <dir> --channels <list> [--tau|--dark|--laser]"
	exit 1
fi

# === VALIDAZIONE INPUT ===
if ! [[ "$N_RUNS" =~ ^[0-9]+$ ]]; then
	echo "❌ --n_runs deve essere un numero intero positivo"
	exit 1
fi
if ! [[ "$VBIAS" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
	echo "❌ --vbias deve essere un numero (es. 10 o 10.5)"
	exit 1
fi
if ! [[ "$SAMPLING" =~ ^[0-9]+$ ]]; then
	echo "❌ --sampling deve essere un numero intero positivo"
	exit 1
fi

# === NORMALIZZA PERCORSI ===
FOLDER="${FOLDER%/}"  # rimuove eventuale slash finale
TARGET_DIR="${FOLDER}/vbias_${VBIAS}"
echo "📁 Output folder: $TARGET_DIR"
mkdir -p "$TARGET_DIR/logs" || { echo "❌ Errore nella creazione di $TARGET_DIR/logs"; exit 1; }

# === TRAP PER CLEANUP ===
cleanup() {
	echo "❌ Script interrotto, cleanup in corso..."
	# Aggiungi qui comandi di cleanup (es. spegnere bias se attivo)
}
trap cleanup EXIT

# === LOOP DEI RUN (run_dgz sequenziali, analisi parallela) ===
all_analysis_pids=()

for ((i=0; i<N_RUNS; i++)); do
    FILENAME="run-$i"
    echo -e "\e[36m>>> Starting acquisition $i with filename: $FILENAME\e[0m"

    SECONDS=0

    # === ACQUISIZIONE DGZ in foreground ===
    echo -e "\e[1m**Lancio run_dgz.py per run $i**\e[0m"
    python run_dgz.py \
        --vbias "$VBIAS" \
        --filter_ADC 0 \
        --filename "$FILENAME" \
        --sampling "$SAMPLING" \
        --trg NIM \
        --channel 0 "${CHANNELS[@]}" \
        --folder "$FOLDER" \
        --min_events "$N_WAVES" || {
        echo -e "\e[31m❌ Errore in run_dgz.py per run $i, esco\e[0m"
        exit 1
    }

    # === LOOP SU CANALI (analisi parallela) ===
    for CH in "${CHANNELS[@]}"; do
        LOGFILE="${TARGET_DIR}/logs/run-${i}_ch${CH}.log"
        touch "$LOGFILE" || { echo "❌ Errore nella creazione di $LOGFILE"; exit 1; }
        FULL_FILENAME="${FOLDER}/${FILENAME}_ch${CH}.npz"
        if [[ ! -f "$FULL_FILENAME" ]]; then
            echo -e "\e[33m⚠️ File $FULL_FILENAME non trovato, salto canale $CH\e[0m"
            continue
        fi

        args=("--filename" "$FULL_FILENAME" "--vbias" "$VBIAS" "--sampling" "$SAMPLING" "--folder" "$FOLDER" "--channel" "$CH")
        $DO_TAU   && args+=("--tau")
        $DO_DARK  && args+=("--dark")
        $DO_LASER && args+=("--laser")

        # Analisi in background
        bash -c "./analyze_channel.sh $(printf "%q " "${args[@]}")" > "$LOGFILE" 2>&1 &
        all_analysis_pids+=($!)
    done

    duration=$SECONDS
    echo -e "\e[35m>>> Run $i completed in $((duration / 60))m $((duration % 60))s\e[0m"
done

# === WAIT GENERALE PER TUTTE LE ANALISI CANALI ===
echo "⏳ Waiting for all channel analyses to complete..."
for pid in "${all_analysis_pids[@]}"; do
    wait "$pid"
done
unset all_analysis_pids
echo -e "\e[32m✅ Tutte le analisi canali completate\e[0m"




# # === SPOSTA FILE NPZ ===
echo "📦 Spostando i file .npz..."
for file in "${FOLDER}"/run-*_*ch*.npz; do
	[[ -f "$file" ]] && mv "$file" "$TARGET_DIR/" || echo "⚠️ Nessun file .npz trovato in ${FOLDER}"
done
echo "✅ Tutti i file .npz spostati in $TARGET_DIR"




# === MERGE DEI FILE (solo se --dark) ===
MERGE_OUTPUT=""  # evita "unbound variable" se non in modalità dark

if $DO_DARK; then
    echo "🧩 Modalità dark: eseguo il merge dei file..."

    MERGE_PATTERN=${FOLDER}/vbias_${VBIAS}_run-*ch*.dark_transition.txt
    MERGE_OUTPUT=${FOLDER}/merged.vbias_${VBIAS}.dark_transitions.txt

    # Fai espandere il glob qui
    files=( $MERGE_PATTERN )
    echo "🔍 File trovati per il merge (pattern: $MERGE_PATTERN):"
    for f in "${files[@]}"; do
        echo "   • $f"
    done
    echo "----------------------------------------"

    if [[ ${#files[@]} -eq 0 ]]; then
        echo -e "\e[33m⚠️ Nessun file trovato per il merge con pattern: $MERGE_PATTERN\e[0m"
    else
        echo -e "\e[1m**Lancio merge.py con ${#files[@]} file:** ${files[*]}\e[0m"

        python ../analysis/script/code/merge.py \
            --merge ${files[*]} \
            --mode sum \
            --columns 2 \
            --sum 2 \
            --output "$MERGE_OUTPUT" || {
            echo -e "\e[31m❌ Errore nel merge dei file dark\e[0m"
            exit 1
        }

        echo -e "\e[32m✅ Merge completato con successo: $MERGE_OUTPUT\e[0m"
    fi
fi

# === ANALISI ROOT solo se il merge è stato eseguito (dark mode) ===
if [[ -n "$MERGE_OUTPUT" && -f "$MERGE_OUTPUT" ]]; then
    echo "📊 Lancio analisi ROOT sul file fuso: $MERGE_OUTPUT"

    ROOT_OUTPUT="${FOLDER}/analysis.vbias_${VBIAS}.root"
    ROOT_MACRO_ABS="/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/root-analysis/rooter_dcr.cpp"

    # Se il file ROOT non esiste, crealo vuoto
    if [[ ! -f "$ROOT_OUTPUT" ]]; then
        echo "🟡 File ROOT non esiste ancora, creo un file vuoto per UPDATE..."
        root -l -b -q <<EOF
TFile f("${ROOT_OUTPUT}", "RECREATE");
f.Close();
.q
EOF
    fi

    # Nome locale della macro (senza VBIAS)
    LOCAL_MACRO="./rooter_dcr_tmp.C"

    echo "🧩 Creo macro temporanea ROOT: $LOCAL_MACRO"
    cat > "$LOCAL_MACRO" <<EOF
#include "${ROOT_MACRO_ABS}"

void rooter_dcr_tmp() {
    std::vector<std::tuple<std::string, std::string, std::string>> inputFiles = {
        {"${MERGE_OUTPUT}", "${VBIAS}", "data"}
    };
    rooter_dcr("${ROOT_OUTPUT}", inputFiles, "Sensor_${VBIAS}");
}
EOF

    echo "📜 Contenuto macro temporanea:"
    cat "$LOCAL_MACRO"
    echo "---------------------------------------"

    # Esegui la macro compilata con ROOT (funzione auto-call)
    echo "🚀 Eseguo: root -l -b -q '$LOCAL_MACRO'"
    root -l -b -q "$LOCAL_MACRO"

    RET=$?
    if [[ $RET -eq 0 ]]; then
        echo -e "\e[32m✅ Analisi completata: ${ROOT_OUTPUT}\e[0m"
        rm -f "$LOCAL_MACRO"
    else
        echo -e "\e[31m❌ Errore ROOT (codice $RET), macro non rimossa per debug: $LOCAL_MACRO\e[0m"
    fi

else
    echo -e "\e[33mℹ️ Nessun file di merge da analizzare con ROOT (no modalità dark)\e[0m"
fi






# === ANALISI LASER (solo se richiesto) ===
if $DO_LASER; then
    echo "💡 Lancio analisi LASER con ROOT..."

    LASER_MACRO="./laser_analysis_tmp.C"
    ROOT_INPUT="${FOLDER}/rooted_${FILENAME}.root"
    ROOT_OUTPUT="${TARGET_DIR}/laser_analysis.vbias_${VBIAS}.root"
    ROOT_CPP="/home/eic/RICCARDO/caen-DT5742/caendgz-sipmanalysis/analysis/root-analysis/timing.cpp"

    echo "🧩 Creo macro ROOT temporanea: $LASER_MACRO"
    cat > "$LASER_MACRO" <<EOF
#include "${ROOT_CPP}"

void laser_analysis_tmp() {
    timing("${ROOT_INPUT}", "${ROOT_OUTPUT}");
}
EOF

    echo "📜 Contenuto macro:"
    cat "$LASER_MACRO"
    echo "---------------------------------------"

    echo "🚀 Eseguo macro ROOT per analisi laser: root -l -b -q $LASER_MACRO"
    root -l -b -q "$LASER_MACRO"

    RET=$?
    if [[ $RET -eq 0 ]]; then
        echo -e "\e[32m✅ Analisi LASER completata: ${ROOT_OUTPUT}\e[0m"
        rm -f "$LASER_MACRO"
    else
        echo -e "\e[31m❌ Errore ROOT nella macro LASER (codice $RET), controlla $LASER_MACRO\e[0m"
    fi
fi



# === RIMUOVI TRAP ===
trap - EXIT
