# Example of deuteron fit usage. Plots an example event from Run Group C.

from deuteron_fit import DFits
import matplotlib.pyplot as plt
import json

with open("example_data/example_data.json","r") as event:
    for line in event:
        json_dict = json.loads(line.rstrip('\n|\r'))
        signal = json_dict['fitsub']
        freq_list = json_dict['freq_list']

params = { 'A': 0.03,
      'G': -0.00003,
      'r': 1.2,   # r greater than 1 is positive polarization
      'wQ': 0.027,
      'wL': 32.69,
      'eta': -0.02,
      'xi': -0.001}

fit_results = DFits(freq_list, signal, params)

r = fit_results.result.params['r'].value
pol = (r * r - 1) / (r * r + r + 1)

fit = fit_results.result.best_fit

plt.plot(freq_list, signal)
plt.plot(freq_list, fit, "-r")
plt.grid()
plt.title(f'NMR Signal, Deuteron Polarization {pol*100:.2f}%')
plt.show()