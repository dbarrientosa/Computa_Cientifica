import numpy as np


class TimeSeries:
    """ Init signature: TimeSeries(self, time, magnitude)
        Type:           type 
        Objeto serie de tiempo, debe ser inicializado con dos numpy arrays del mismo tamaño. Guarda la curva de luz estandarizada,
        magnitud promedio 0 y desviación estándar 1.    
    """
    
    def __init__(self, time, magnitude):
        
        self.mag = (magnitude - magnitude.mean()) / magnitude.std()
        self.time = time
       
    def subsample(self, ti, l):
    """ Función subsample(self, t_i, l):
    Extracción de un subsample de la curva de luz. Se recibe un tiempo inicial y un intervalo de tiempo a samplear. El intervalo
    de tiempo entregado comienza desde 0.

    input:  t_i, tiempo inicial a samplear, debe estar contenido en el intervalo de la curva.
            l, intervalo de tiempo a samplear, la curva debe tener al menos largo t_i + l.
    output: subsample, Objeto TimeSeries(), el subsample de la curva.
    """
        #Identificar indices de tiempo inicial y tiempo final
        i_ti = np.argmin(np.abs(self.time - ti))
        i_tf = np.argmin(np.abs(self.time - ti - l))
        
        #Definir nuevos intervalos time (desde 0) y mag
        time_new = self.time[i_ti:i_tf + 1] - self.time[i_ti]
        mag_new = self.mag[i_ti:i_tf + 1]
        
        
        return TimeSeries(time_new, mag_new)

    
        
    
