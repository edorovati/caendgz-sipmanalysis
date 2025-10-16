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

# === LOOP DEI RUN ===
for ((i=0; i<N_RUNS; i++)); do
	FILENAME="run-$i"
	echo -e "\e[36m>>> Starting acquisition $i with filename: $FILENAME\e[0m"

	SECONDS=0

	# === ACQUISIZIONE DGZ ===
	echo -e "\e[1m**Lancio run_dgz.py con argomenti:** --vbias \"$VBIAS\" --filter_ADC 0 --filename \"$FILENAME\" --sampling \"$SAMPLING\" --trg NIM --channel 0 ${CHANNELS[*]} --folder \"$FOLDER\" --min_events \"$N_WAVES\"\e[0m"
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

	# === LOOP SU CANALI ===
	pids=()
	for CH in "${CHANNELS[@]}"; do
		LOGFILE="${TARGET_DIR}/logs/run-${i}_ch${CH}.log"
		touch "$LOGFILE" || { echo "❌ Errore nella creazione di $LOGFILE"; exit 1; }
		echo -e "\e[32m>>> [RUN $i | ch$CH] Starting analysis (log: $LOGFILE)\e[0m"

		FULL_FILENAME="${FOLDER}/${FILENAME}_ch${CH}.npz"
		# Verifica che il file esista
		if [[ ! -f "$FULL_FILENAME" ]]; then
			echo -e "\e[31m❌ File $FULL_FILENAME non trovato, salto analisi per canale $CH\e[0m"
			continue
		fi
		args=("--filename" "$FULL_FILENAME" "--vbias" "$VBIAS" "--sampling" "$SAMPLING" "--folder" "$FOLDER" "--channel" "$CH")
		$DO_TAU   && args+=("--tau")
		$DO_DARK  && args+=("--dark")
		$DO_LASER && args+=("--laser")

		# Stampa argomenti per analyze_channel.sh
		args_str=$(printf "%q " "${args[@]}")
		echo -e "\e[1m**Lancio analyze_channel.sh con argomenti:** $args_str\e[0m"

		# Esegui analyze_channel.sh con quoting migliorato
		bash -c "./analyze_channel.sh $(printf "%q " "${args[@]}")" > "$LOGFILE" 2>&1 &

		pids+=($!)
	done

	# === WAIT e controllo errori ===
	echo "⏳ Waiting for analyses of run $i..."
	for idx in "${!pids[@]}"; do
		pid=${pids[$idx]}
		CH=${CHANNELS[$idx]}
		LOGFILE="${TARGET_DIR}/logs/run-${i}_ch${CH}.log"
		if ! wait "$pid"; then
			echo -e "\e[31m❌ ERROR detected in run $i, channel $CH. Check log: $LOGFILE\e[0m"
			grep '❌' "$LOGFILE" | while read -r line; do
				echo -e "\e[31m  $line\e[0m"
			done
		else
			echo -e "\e[32m✅ Analysis for run $i, channel $CH completed successfully\e[0m"
		fi
	done
	unset pids

	duration=$SECONDS
	echo -e "\e[35m>>> Run $i completed in $((duration / 60))m $((duration % 60))s\e[0m"
done

# # === SPOSTA FILE NPZ ===
echo "📦 Spostando i file .npz..."
for file in "${FOLDER}"/run-*_*ch*.npz; do
	[[ -f "$file" ]] && mv "$file" "$TARGET_DIR/" || echo "⚠️ Nessun file .npz trovato in ${FOLDER}"
done
echo "✅ Tutti i file .npz spostati in $TARGET_DIR"

# === RIMUOVI TRAP ===
trap - EXIT