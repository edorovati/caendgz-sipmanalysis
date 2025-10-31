#!/bin/bash
# Script per controllare l'alimentatore Aim-TTi PLH120-P tramite TCP/IP

host="aimtti-plh120p-00"
port=9221

usage() {
    echo "Usage: $0 --vbias <value> [--on | --off]"
    echo "  --vbias   Valore della tensione di bias (in Volt)"
    echo "  --on      Accende l'uscita (OP1 1)"
    echo "  --off     Spegne l'uscita (OP1 0)"
    exit 1
}

# Parse argomenti
if [[ $# -lt 2 ]]; then
    usage
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vbias)
            VBIAS="$2"
            shift 2
            ;;
        --on)
            ACTION="on"
            shift
            ;;
        --off)
            ACTION="off"
            shift
            ;;
        *)
            usage
            ;;
    esac
done

# Controllo argomenti obbligatori
if [[ -z "$VBIAS" ]]; then
    echo "Errore: devi specificare --vbias <value>"
    usage
fi

if [[ -z "$ACTION" ]]; then
    echo "Errore: devi specificare --on oppure --off"
    usage
fi

# Esecuzione comando
echo "Setting V1 to ${VBIAS} V..."
echo "V1 $VBIAS" | nc -w1 -W1 "$host" "$port"

if [[ "$ACTION" == "on" ]]; then
    echo "Switching output ON..."
    echo "OP1 1" | nc -w1 -W1 "$host" "$port"
    sleep 10
    echo "Output is ON and V1 is set to ${VBIAS} V."
else
    echo "Switching output OFF..."
    echo "OP1 0" | nc -w1 -W1 "$host" "$port"
    echo "Output is OFF."
fi
