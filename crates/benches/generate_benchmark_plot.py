import json, os, glob
import matplotlib.pyplot as plt

# Read in files
result = (y for x in os.walk("target/criterion") for y in glob.glob(os.path.join(x[0], 'estimates.json')))

rustfft = {}
fftw = {}
monarch = {}

for file in result:
    parts = file.split('/')
    test = parts[2].split('-')
    length = int(test[1])
    with open(file, 'r') as f:
        results = json.load(f)
    mean = results['mean']['point_estimate']

    if 'monarch' in test[0]:
        monarch[length] = mean
    elif 'fftw' in test[0]:
        fftw[length] = mean
    elif 'rustfft' in test[0]:
        rustfft[length] = mean

x = list(rustfft.keys())
x.sort()

rustfft = [rustfft[i] for i in x]
fftw = [fftw[i] for i in x]
monarch = [monarch[i] for i in x]

plt.plot(x, rustfft, label='rustfft')
plt.plot(x, fftw, label='fftw')
plt.plot(x, monarch, label='monarch')
plt.title('FFT Times')
plt.ylabel('Time (ns)')
plt.xlabel('FFT Size')
plt.legend(loc='upper left')
plt.show()
