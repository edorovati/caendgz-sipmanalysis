#ifndef TIMESTUFF_H
#define TIMESTUFF_H
#include <vector>
class TimeStuff {
public:
    // Calcola il crossing con interpolazione lineare
    static double computeCrossingTime(const Double_t* time, const Double_t* amp, int size, double threshold) {
        for (int i = 1; i < size; ++i) {
            if (amp[i - 1] < threshold && amp[i] >= threshold) {
                double x0 = time[i - 1], x1 = time[i];
                double y0 = amp[i - 1], y1 = amp[i];
                return x0 + (threshold - y0) * (x1 - x0) / (y1 - y0);
            }
        }
        return -1; // Crossing non trovato
    }
    static double computeCrossingTimeAboveDelta(const Double_t* time, const Double_t* amp, int size, double threshold, double crossing_time_laser, double min_delta = 30.0) {
        for (int i = 1; i < size; ++i) {
            if (amp[i - 1] < threshold && amp[i] >= threshold) {
                double x0 = time[i - 1], x1 = time[i];
                double y0 = amp[i - 1], y1 = amp[i];
                double crossing_time = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0);
                if ((crossing_time - crossing_time_laser) > min_delta) {
                    return crossing_time;
                }
            }
        }
        return -1; // Crossing non trovato con delta > min_delta
    }


    
    // Funzione per calcolare tutti i crossing per una soglia
    static std::vector<double> computeAllCrossingTimes(const Double_t* time, const Double_t* amp, int size, double threshold) {
        std::vector<double> crossings;
        for (int i = 1; i < size; ++i) {
            if (amp[i - 1] < threshold && amp[i] >= threshold) {
                double x0 = time[i - 1], x1 = time[i];
                double y0 = amp[i - 1], y1 = amp[i];
                double cross_time = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0);
                crossings.push_back(cross_time);
            }
        }
        return crossings;
    }

    
    static int findMaxAfterCrossing(const std::vector<double>& amps,
                                int crossing_index,
                                double threshold)
    {
        if (crossing_index < 0 || crossing_index >= static_cast<int>(amps.size()))
            return -1; // indice non valido
        
        int max_idx = crossing_index;
        double max_val = amps[crossing_index];
        double hysteresis = threshold * 0.5;  // isteresi = 50% della soglia
        
        for (size_t i = crossing_index + 1; i < amps.size(); ++i) {
            if (amps[i] > max_val) {
                max_val = amps[i];
                max_idx = static_cast<int>(i);
            }
            // Fermati solo quando scende di almeno hysteresis rispetto al massimo
            if (amps[i] < max_val - hysteresis) {
                break;
            }
        }
        
        // Se il massimo non supera la soglia richiesta, scartalo
        if (max_val < threshold) {
            return -1;
        }
        
        return max_idx;
    }

    /*struct CrossingResult {
        int max_index;        // indice del massimo dopo il crossing
        double crossing_time; // tempo del crossing (interpolato)
    };

    // Funzione "suina": trova il crossing, il massimo locale e scarta waveform sporche
    static CrossingResult analyzeWaveform(const Double_t* time,
                                          const Double_t* amp,
                                          int size,
                                          double threshold)
    {
        CrossingResult result{-1, -1.0};

        // 1) Trova crossing con interpolazione lineare
        int crossing_index = -1;
        for (int i = 1; i < size; ++i) {
            if (amp[i - 1] < threshold && amp[i] >= threshold) {
                double x0 = time[i - 1], x1 = time[i];
                double y0 = amp[i - 1], y1 = amp[i];
                result.crossing_time = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0);
                crossing_index = i;
                break;
            }
        }

        if (crossing_index == -1) {
            return result; // nessun crossing trovato
        }

        // 2) Cerca il massimo locale: fermati quando trovi due punti consecutivi decrescenti
        int max_idx = crossing_index;
        double max_val = amp[crossing_index];

        for (int i = crossing_index + 1; i < size - 2; ++i) {
            if (amp[i] > max_val) {
                max_val = amp[i];
                max_idx = i;
            }
            // due punti consecutivi in discesa → fermo la ricerca del massimo
            if (amp[i+1] < amp[i] && amp[i+2] < amp[i+1]) {
                break;
            }
        }

        // 3) Controlla il resto della WF: se esiste un valore più grande → scarta (tempo invalido)
        for (int i = max_idx + 1; i < size; ++i) {
            if (amp[i] > max_val) {
                return CrossingResult{-1, -1.0}; // waveform sporca
            }
        }

        // 4) Se il massimo non supera la soglia → scarta
        if (max_val < threshold) {
            return CrossingResult{-1, -1.0};
        }

        result.max_index = max_idx;
        return result;
    }*/
    
    
 // Funzione combinata
    /*static CrossingResult analyzeWaveform(const Double_t* time,
                                       const Double_t* amp,
                                       int size,
                                       double threshold)
    {
        CrossingResult result{-1, -1.0};
        
        // 1) Trova crossing con interpolazione lineare
        int crossing_index = -1;
        for (int i = 1; i < size; ++i) {
            if (amp[i - 1] < threshold && amp[i] >= threshold) {
                double x0 = time[i - 1], x1 = time[i];
                double y0 = amp[i - 1], y1 = amp[i];
                result.crossing_time = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0);
                crossing_index = i;
                break;
            }
        }
        
        if (crossing_index == -1) {
            return result; // nessun crossing trovato
        }
        
        // 2) Trova massimo dopo crossing con isteresi = 50% della soglia
        int max_idx = crossing_index;
        double max_val = amp[crossing_index];
        double hysteresis = threshold * 0.5;
        
        for (int i = crossing_index + 1; i < size; ++i) {
            if (amp[i] > max_val) {
                max_val = amp[i];
                max_idx = i;
            }
            if (amp[i] < max_val - hysteresis) {
                break;
            }
        }
        
        // Se non supera la soglia → scarta
        if (max_val < threshold) {
            return result;
        }
        
        result.max_index = max_idx;
        return result;
    }*/
    
    /*static CrossingResult analyzeWaveform(const Double_t* time,
                                          const Double_t* amp,
                                          int size,
                                          double threshold)
    {
        CrossingResult result{-1, -1.0};
        
        // 1) Trova crossing con interpolazione lineare
        int crossing_index = -1;
        for (int i = 1; i < size; ++i) {
            if (amp[i - 1] < threshold && amp[i] >= threshold) {
                double x0 = time[i - 1], x1 = time[i];
                double y0 = amp[i - 1], y1 = amp[i];
                result.crossing_time = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0);
                crossing_index = i;
                break;
            }
        }
        
        if (crossing_index == -1) {
            return result; // nessun crossing trovato
        }
        
        // 2) Trova massimo dopo crossing
        int max_idx = crossing_index;
        double max_val = amp[crossing_index];
        
        int dec_count = 0;               // contatore punti decrescenti consecutivi
        const int DEC_LIMIT = 30;        // primi 15 punti decrescenti
        
        for (int i = crossing_index + 1; i < size; ++i) {
            if (amp[i] > max_val) {
                // Nuovo massimo → reset contatore
                max_val = amp[i];
                max_idx = i;
                dec_count = 0;
            } else {
                // Punto decrescente
                dec_count++;
            }

            // Se abbiamo 15 punti decrescenti consecutivi → stop
            if (dec_count >= DEC_LIMIT) {
                break;
            }
        }
        
        // Se non supera la soglia → scarta
        if (max_val < threshold) {
            return result;
        }
        
        result.max_index = max_idx;
        return result;
    }*/
    
    /*static CrossingResult analyzeWaveform(const Double_t* time,
                                          const Double_t* amp,
                                          int size,
                                          double threshold)
    {
        CrossingResult result{-1, -1.0};
        
        // 1) Trova crossing con interpolazione lineare
        int crossing_index = -1;
        for (int i = 1; i < size; ++i) {
            if (amp[i - 1] < threshold && amp[i] >= threshold) {
                double x0 = time[i - 1], x1 = time[i];
                double y0 = amp[i - 1], y1 = amp[i];
                result.crossing_time = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0);
                crossing_index = i;
                break;
            }
        }
        
        if (crossing_index == -1) {
            return result; // nessun crossing trovato
        }
        
        // 2) Trova massimo assoluto dopo crossing (fino alla fine del waveform)
        int max_idx = crossing_index;
        double max_val = amp[crossing_index];
        
        for (int i = crossing_index + 1; i < size; ++i) {
            if (amp[i] > max_val) {
                max_val = amp[i];
                max_idx = i;
            }
        }
        
        // Se non supera la soglia → scarta
        if (max_val < threshold) {
            return result;
        }
        
        result.max_index = max_idx;
        return result;
    }*/
    struct CrossingResult {
        int max_index;          // indice del massimo dopo il crossing
        double crossing_time;   // tempo del crossing (interpolato)
    };

    // Nuova versione di analyzeWaveform
    static CrossingResult analyzeWaveform(const Double_t* time,
                                          const Double_t* amp,
                                          int size,
                                          double threshold)
    {
        CrossingResult result{-1, -300.0};  // inizializzo crossing_time a -300 ns

        // trova crossing lineare con controllo di salita
        int crossing_index = -1;
        for (int i = 1; i < size - 3; ++i) {  // -3 perché guardiamo i+3
            if (amp[i] < threshold && amp[i+1] >= threshold) {
                // controllo che il fronte sia in crescita
                if (amp[i-1] < amp[i] && amp[i+1] < amp[i+2] && amp[i+2] < amp[i+3]) {
                    double x0 = time[i], x1 = time[i+1];
                    double y0 = amp[i], y1 = amp[i+1];
                    result.crossing_time = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0);
                    crossing_index = i;
                    break;
                }
            }
        }

        if (crossing_index == -1) return result;  // nessun crossing valido

        // cerca massimo entro 3 ns dal crossing
        int max_idx = crossing_index;
        double max_val = amp[crossing_index];
        double t_limit = result.crossing_time + 3.0;

        for (int i = crossing_index + 1; i < size; ++i) {
            if (time[i] > t_limit) break;
            if (amp[i] > max_val) {
                max_val = amp[i];
                max_idx = i;
            }
        }

        result.max_index = max_idx;
        return result;
    }
    
    // CFD: Constant Fraction Discriminator
    // Restituisce crossing_time al "fraction" del massimo (es. 0.5 = 50%)
    // max_idx deve essere già stato trovato da analyzeWaveform
    static CrossingResult computeCFDCrossing(const Double_t* time,
                                             const Double_t* amp,
                                             int size,
                                             int max_idx,
                                             double fraction = 0.5)
    {
        CrossingResult result{-1, -300.0};  // init crossing_time a -300 ns

        if (size < 2 || max_idx <= 0 || max_idx >= size)
            return result;

        double max_amp = amp[max_idx];
        double threshold = fraction * max_amp;

        int crossing_index = -1;

        // Scansiona in avanti per trovare il crossing della soglia
        for (int i = 1; i < max_idx; ++i) {
            if (amp[i] < threshold && amp[i+1] >= threshold) {
                // Controllo che il fronte sia in salita
                if (amp[i-1] < amp[i] && amp[i+1] < amp[i+2] && amp[i+2] < amp[i+3]) {
                    double x0 = time[i], x1 = time[i+1];
                    double y0 = amp[i], y1 = amp[i+1];
                    result.crossing_time = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0);
                    crossing_index = i;
                    break;
                }
            }
        }

        if (crossing_index == -1) return result;  // nessun crossing valido

        // cerca massimo locale vicino al crossing (entro 3 ns)
        int local_max_idx = crossing_index;
        double max_val = amp[crossing_index];
        double t_limit = result.crossing_time + 3.0;

        for (int i = crossing_index + 1; i < size; ++i) {
            if (time[i] > t_limit) break;
            if (amp[i] > max_val) {
                max_val = amp[i];
                local_max_idx = i;
            }
        }

        // se il massimo trovato non coincide con max_idx, continua a cercare fino a max_idx
        if (local_max_idx != max_idx) {
            for (int i = local_max_idx + 1; i <= max_idx; ++i) {
                if (amp[i] > max_val) {
                    max_val = amp[i];
                    local_max_idx = i;
                }
            }
        }

        result.max_index = local_max_idx;
        return result;
    }

    static std::pair<double, double> computeLinearFitWithSlope(const Double_t* time,
                                                               const Double_t* amp,
                                                               int size,
                                                               int max_idx)
    {
        if (size < 2 || max_idx <= 0 || max_idx >= size)
            return {-1.0, 0.0};

        double max_amp = amp[max_idx];
        double thr80 = 0.8 * max_amp;

        // --- crossing 80% ---
        double t80 = -1.0;
        for (int i = 1; i <= max_idx; ++i) {
            if (amp[i - 1] < thr80 && amp[i] >= thr80) {
                double x0 = time[i - 1], x1 = time[i];
                double y0 = amp[i - 1], y1 = amp[i];
                t80 = x0 + (thr80 - y0) * (x1 - x0) / (y1 - y0);
                break;
            }
        }
        if (t80 < 0) return {-1.0, 0.0};

        // --- crossing 20% (a ritroso) ---
        double thr20 = 0.2 * max_amp;
        int idx20 = -1;
        for (int i = 1; i <= max_idx; ++i) {
            int j = max_idx - i;
            if (j <= 0) break;
            if (amp[j - 1] < thr20 && amp[j] >= thr20) {
                idx20 = j;
                break;
            }
        }
        if (idx20 < 0) return {-1.0, 0.0};

        double x0 = time[idx20 - 1], x1 = time[idx20];
        double y0 = amp[idx20 - 1], y1 = amp[idx20];
        double t20 = x0 + (thr20 - y0) * (x1 - x0) / (y1 - y0);

        // --- retta sul fronte di salita ---
        double slope = (thr80 - thr20) / (t80 - t20);
        double intercept = thr20 - slope * t20;

        // intercetta sull’asse y=0
        double t_cross = -intercept / slope;

        return {t_cross, slope};
    }

    
    // Versione dopo 35 ns
    static CrossingResult analyzeWaveformAfter35ns(const Double_t* time,
                                                   const Double_t* amp,
                                                   int size,
                                                   double threshold)
    {
        CrossingResult result{-1, -300.0};  // crossing_time inizializzato a -300 ns

        int crossing_index = -1;
        for (int i = 1; i < size-1; ++i) {
            if (time[i] < 38.0) continue;
            if (time[i] > 42.0) break;

            if (amp[i] < threshold && amp[i+1] >= threshold) {
                double x0 = time[i], x1 = time[i+1];
                double y0 = amp[i], y1 = amp[i+1];
                result.crossing_time = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0);
                crossing_index = i;
                break;
            }
        }

        if (crossing_index == -1) return result;

        // cerca massimo fino a quando i 3 punti successivi sono decrescenti
        int max_idx = crossing_index;
        double max_val = amp[crossing_index];

        for (int i = crossing_index + 1; i < size - 3; ++i) {
            if (amp[i] > max_val) {
                max_val = amp[i];
                max_idx = i;
            }

            if (amp[i] > amp[i+1] && amp[i+1] > amp[i+2] && amp[i+2] > amp[i+3]) {
                break;
            }
        }

        result.max_index = max_idx;
        return result;
    }



};

#endif // TIMESTUFF_H
