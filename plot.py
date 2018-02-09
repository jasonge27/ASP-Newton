import matplotlib.pyplot as plt
import numpy as np

plt.style.use('ggplot')

dataset_name = "simic"

f1 = open(dataset_name + "_out.txt", 'r')
lines = f1.readlines()
f1.close()

lines = [line.strip() for line in lines]

optimal_obj = 0.14

ASPNewton_time = np.array([float(x) for x in lines[1].split(',')])
ASPNewton_obj = np.array([float(x) for x in lines[2].split(',')])

glmnet_time = np.array([float(x) for x in lines[4].split(',')])
glmnet_obj = np.array([float(x) for x in lines[5].split(',')])

gcdnet_time = np.array([float(x) for x in lines[7].split(',')])
gcdnet_obj = np.array([float(x) for x in lines[8].split(',')])

fista_time = np.array([float(x) for x in lines[10].split(',')])
fista_obj = np.array([float(x) for x in lines[11].split(',')])

f1 = open(dataset_name + "_picasso.txt", 'r')
lines = f1.readlines()
f1.close()
lines = [line.strip() for line in lines]
picasso_time = np.array([float(x) for x in lines[1].split(',')])
picasso_obj = np.array([float(x) for x in lines[2].split(',')])

plt.plot(ASPNewton_time, ASPNewton_obj, marker='o', color='#ff66cc')
plt.plot(picasso_time, picasso_obj, marker='o', color='#4da6ff')
plt.plot(glmnet_time, glmnet_obj, marker='o', color='#ffff14')
plt.plot(gcdnet_time, gcdnet_obj, marker='o', color='#029386')
plt.plot(fista_time, fista_obj, marker='o', color='#001146')

plt.legend(
    ['Ours', 'picasso', 'proximal Newton', 'naive CD', 'accelerated CD'],
    loc=4,
    fontsize=10)

plt.xlabel('CPU Time (s)', fontsize=12)
plt.ylabel('Logarithm Objective Gap', fontsize=12)
plt.title('Convergence Rate on Dataset %s' % (dataset_name), fontsize=12)
plt.savefig('CPUTime_%s.eps' % (dataset_name), format='eps')
#(Warmstart from A Pathwise Solution)
# for lambda =0.02 from
