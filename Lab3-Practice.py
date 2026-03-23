
import pandas
import matplotlib.pyplot as plt
from pandas.tools.plotting import scatter_matrix
import numpy as np




def computeCost(X, y, theta):
	"""
	   computes the cost of using theta as the parameter for linear
	   regression to fit the data points in X and y
	"""

	m = y.size
	J = 0

	h = np.dot(X,theta)
	sq_error = np.sum(np.square(h - y))
	J =  (sq_error) / (2 * m)
# =========================================================================

	return J

def gradient_descent(X, y, theta, alpha, num_iters):
	'''
	Performs gradient descent to learn theta
	by taking num_items gradient steps with learning
	rate alpha and update theata0, theata1
	'''
	for i in num :
		h=theta[0]+theta[1]*X
		r=y.size()
		j=((1/r))*((h-y)**2)
		d0=j
		d1=j*X
		theta[0]=theta[0]-alpha*d0
		theta[1] = theta[1] - alpha * d1
	return theta


#Load the dataset
data = pandas.read_csv('data.csv')

scatter_matrix(data[['population','profit']])
plt.show()

X = data['population']
y = data['profit']

#number of training samples
m = y.size
#X = np.vstack(zip(np.ones(m),data['population']))

#Initialize theta parameters
theta0 = 0
theta1 =0
#Some gradient descent settings
iterations = 1500
alpha = 0.01
#compute and display initial cost
theta = gradient_descent(X, y, [theta0,theta1], alpha, iterations) #To be completed by students
print "theta",theta

#Predict values for population sizes of 3.5 and 7.0
#Students write prediction code

h = np.dot([1,3.5], theta)
print(h)
