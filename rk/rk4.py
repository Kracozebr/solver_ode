import numpy as np
import logging

# Configure logging
logging.basicConfig(filename='rk4.log', level=logging.INFO, 
                    format='%(asctime)s - %(message)s')

class RungeKutta4:
    def __init__(self, system, initial_conditions, t0, tf, h):
        self.system = system
        self.y0 = np.array(initial_conditions)
        self.t0 = t0
        self.tf = tf
        self.h = h

    def rk4_step(self, t, y):
        k1 = self.h * self.system(t, y)
        k2 = self.h * self.system(t + 0.5 * self.h, y + 0.5 * k1)
        k3 = self.h * self.system(t + 0.5 * self.h, y + 0.5 * k2)
        k4 = self.h * self.system(t + self.h, y + k3)
        step_result = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        logging.info(f"t = {t}, y = {y}, k1 = {k1}, k2 = {k2}, k3 = {k3}, k4 = {k4}, step_result = {step_result}")
        return step_result

    def solve(self):
        t_values = np.arange(self.t0, self.tf, self.h)
        y_values = np.zeros((len(t_values), len(self.y0)))
        y_values[0] = self.y0
        for i in range(1, len(t_values)):
            y_values[i] = self.rk4_step(t_values[i-1], y_values[i-1])
        return t_values, y_values
