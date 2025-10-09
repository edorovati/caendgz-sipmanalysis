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
    
    struct CrossingResult {
        int max_index;      // indice del massimo dopo il crossing
        double crossing_time; // tempo del crossing (interpolato)
    };

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
    
    // Nuova versione di analyzeWaveform
    static CrossingResult analyzeWaveform(const Double_t* time,
                                          const Double_t* amp,
                                          int size,
                                          double threshold)
    {
        CrossingResult result{-1, -1.0};

        // trova crossing lineare
        int crossing_index = -1;
        for (int i = 1; i < size; ++i) {
            if (amp[i-1] < threshold && amp[i] >= threshold) {
                double x0 = time[i-1], x1 = time[i];
                double y0 = amp[i-1], y1 = amp[i];
                result.crossing_time = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0);
                crossing_index = i;
                break;
            }
        }

        if (crossing_index == -1) return result;

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



};

#endif // TIMESTUFF_H
