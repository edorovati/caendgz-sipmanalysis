#!/bin/bash
# Script per controllare Aim-TTi PLH120-P via TCP/IP

host="aimtti-plh120p-00"
port=9221

usage() {
    echo "Usage:"
    echo "  $0 --vbias <V> [--on | --off]"
    echo "  $0 --read"
    echo ""
    echo "Options:"
    echo "  --vbias <V>   imposta la tensione in Volt"
    echo "  --on          accende l'uscita"
    echo "  --off         spegne l'uscita"
    echo "  --read        legge tensione, corrente e stato"
    exit 1
}

# Nessun argomento
if [[ $# -lt 1 ]]; then
    usage
fi

# Parsing argomenti
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
        --read)
            ACTION="read"
            shift
            ;;
        *)
            usage
            ;;
    esac
done

# Lettura stato
if [[ "$ACTION" == "read" ]]; then
    echo "📟 Lettura stato alimentatore..."
    echo -n "Tensione impostata: "
    echo "V1?" | nc -w1 -W1 "$host" "$port"
    echo -n "Tensione erogata:   "
    echo "V1O?" | nc -w1 -W1 "$host" "$port"
    echo -n "Corrente erogata:   "
    echo "I1O?" | nc -w1 -W1 "$host" "$port"
    echo -n "Uscita:             "
    STATE=$(echo "OP1?" | nc -w1 -W1 "$host" "$port")
    if [[ "$STATE" == "1" ]]; then
        echo "ON"
    else
        echo "OFF"
    fi
    exit 0
fi

# Controllo parametri
if [[ -z "$VBIAS" ]]; then
    echo "Errore: devi specificare --vbias <valore>"
    usage
fi
if [[ -z "$ACTION" ]]; then
    echo "Errore: devi specificare --on o --off"
    usage
fi

# Imposta tensione
echo "⚙️  Imposto tensione a ${VBIAS} V..."
echo "V1 $VBIAS" | nc -w1 -W1 "$host" "$port"

# Accendi/spegni
if [[ "$ACTION" == "on" ]]; then
    echo "🔌 Accendo l’uscita..."
    echo "OP1 1" | nc -w1 -W1 "$host" "$port"
    sleep 1
    echo "✅ Uscita ON – V=${VBIAS} V"
elif [[ "$ACTION" == "off" ]]; then
    echo "⏻ Spengo l’uscita..."
    echo "OP1 0" | nc -w1 -W1 "$host" "$port"
    echo "✅ Uscita OFF"
fi
