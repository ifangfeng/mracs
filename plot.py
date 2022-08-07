import numpy as np
import matplotlib.pyplot as plt

x = [0.5, 0.523564, 0.548239, 0.574077, 0.601132, 0.629463, 0.659128, 0.690192, 0.72272, 0.756781, 0.792447, 0.829793, 0.8689, 0.90985, 0.95273, 0.997631, 1.04465, 1.09388, 1.14543, 1.19942, 1.25594, 1.31513, 1.37711, 1.44202, 1.50998, 1.58114, 1.65566, 1.73368, 1.81539, 1.90095, 1.99054, 2.08435, 2.18258, 2.28544, 2.39315, 2.50594, 2.62404, 2.7477, 2.8772, 3.0128, 3.15479, 3.30347, 3.45915, 3.62218, 3.79289, 3.97164, 4.15882, 4.35482, 4.56005, 4.77496, 5, 5.23564, 5.48239, 5.74077, 6.01132, 6.29463, 6.59128, 6.90192, 7.2272, 7.56781, 7.92447, 8.29793, 8.689, 9.0985, 9.5273, 9.97631, 10.4465, 10.9388, 11.4543, 11.9942, 12.5594, 13.1513, 13.7711, 14.4202, 15.0998, 15.8114, 16.5566, 17.3368, 18.1539, 19.0095, 19.9054, 20.8435, 21.8258, 22.8544, 23.9315, 25.0594, 26.2404, 27.477, 28.772, 30.128, 31.5479, 33.0347, 34.5915, 36.2218, 37.9289, 39.7164, 41.5882, 43.5482, 45.6005, 47.7496, 50]
y = [352.919, 243.791, 177.035, 134.114, 105.14, 84.6657, 69.5937, 58.103, 49.0886, 41.856, 35.9514, 31.0664, 26.9824, 23.5392, 20.616, 18.1197, 15.9768, 14.129, 12.5292, 11.139, 9.92671, 8.86621, 7.93572, 7.11696, 6.3946, 5.75564, 5.18909, 4.68556, 4.23706, 3.83673, 3.47868, 3.15782, 2.86977, 2.61071, 2.37734, 2.16676, 1.97647, 1.80426, 1.64818, 1.50654, 1.37784, 1.26075, 1.15409, 1.05684, 0.968063, 0.886944, 0.812752, 0.744835, 0.682609, 0.625552, 0.573198, 0.525124, 0.480955, 0.440351, 0.403006, 0.368644, 0.337018, 0.307905, 0.281102, 0.256428, 0.233717, 0.21282, 0.193601, 0.175934, 0.159705, 0.144808, 0.131146, 0.118628, 0.107169, 0.0966924, 0.0871253, 0.0784008, 0.0704569, 0.0632365, 0.0566865, 0.0507574, 0.045403, 0.0405792, 0.0362443, 0.0323579, 0.028881, 0.0257761, 0.0230067, 0.0205381, 0.0183375, 0.0163742, 0.0146201, 0.01305, 0.0116416, 0.0103753, 0.00923467, 0.00820555, 0.00727608, 0.00643628, 0.00567772, 0.00499317, 0.00437638, 0.00382184, 0.00332461, 0.00288015, 0.00248426]
plt.scatter(x, y)
plt.semilogx(base = 10)
plt.semilogy(base = 10)
plt.title("Millennium-I Galaxy")
plt.xlabel("r Mpc/h")
plt.ylabel(r'$\xi$(r)')
plt.show()