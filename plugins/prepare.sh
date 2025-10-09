#!/bin/bash

echo "🔧 [prepare.sh] Preparazione ambiente..."

# Attiva ambiente virtuale Python
echo "🧪 Attivazione ambiente virtuale..."
source ~/RICCARDO/myenv/bin/activate

# Riavvia server
echo "🛑 Kill vecchi server..."
./kill_servers.sh
sleep 1

echo "🚀 Avvio nuovi server..."
./start_servers.sh

# Aspetta che aimtti-server.py parta
for i in {1..10}; do
    if pgrep -f aimtti-server.py > /dev/null; then
        echo "✅ Processo aimtti-server.py attivo"
        break
    fi
    echo "⏳ Attendo che aimtti-server.py si avvii..."
    sleep 0.5
done

# Cerca la porta TCP su cui sta ascoltando
PORT=$(sudo lsof -i -nP | grep LISTEN | grep aimtti | awk '{print $9}' | sed 's/.*://')

# # Se trovata, attendi che sia attiva
# if [[ -n "$PORT" ]]; then
#     for i in {1..10}; do
#         if nc -z localhost "$PORT"; then
#             echo "✅ Server pronto su porta $PORT"
#             break
#         fi
#         echo "⏳ Attesa porta $PORT..."
#         sleep 0.5
#     done
# else
#     echo "⚠️  Porta del server non trovata. Continuo comunque..."
#     sleep 1
# fi

echo "✅ Servers pronti."