import pandas as pd

# problem 1 get the version of pandas
# print(pd.__version__)

# problem 2 convert to series data in pandas
# stocks = ['PLW', 'CDR', '11B', 'TEN']
# print(pd.Series(data=stocks))

# problem 3 print series of data set with pandas
# stocks = {'PLW': 387.00, 'CDR': 339.5, 'TEN': 349.5, '11B': 391.0}
# print(pd.Series(stocks))

# problem 4 convert it to a list
# stocks = {'PLW': 387.00, 'CDR': 339.5, 'TEN': 349.5, '11B': 391.0}
# quotations = pd.Series(data=stocks).tolist()
# print(quotations)

# problem 5 name column of dataset
# stocks = {'PLW': 387.00, 'CDR': 339.5, 'TEN': 349.5, '11B': 391.0}
# quotations = pd.Series(data=stocks)
# quotations =pd.DataFrame(quotations,columns=['price'])
# print(quotations)

# problem 6 make a range by numpy and print it
# print(pd.Series(data=np.arange(10, 100, 10), index=np.arange(101, 110), dtype='float'))

# problem 7 convert to type int
# series = pd.Series(['001', '002', '003', '004'], list('abcd')).astype(int)
# print(series)
# or
# series = pd.to_numeric(series)
# print(series)

# problem 8
# stocks = {'PLW': 387.00, 'CDR': 339.5, 'TEN': 349.5, '11B': 391.0}
# newStock={'BBT': 25.5, 'F51': 19.2}
# # quotations = pd.Series(stocks).append(pd.Series(newStock)) -> Future Warning instead using concat()
# quotations = pd.concat([pd.Series(stocks),pd.Series(newStock)])
# print(quotations)

# problem 9
# stocks = {
#     'PLW': 387.00,
#     'CDR': 339.5,
#     'TEN': 349.5,
#     '11B': 391.0,
#     'BBT': 25.5,
#     'F51': 19.2
# }
# quotations = pd.Series(data=stocks)
# quotations = pd.DataFrame(quotations).reset_index()
# quotations.columns = ['ticker','price']
# print(quotations)


# problem 10
# data_dict = {
#     'company': ['Amazon', 'Microsoft', 'Facebook'],
#     'price': [2375.00, 178.6, 179.2],
#     'ticker': ['AMZN.US', 'MSFT.US', 'FB.US']
# }
# companies = pd.DataFrame(data=data_dict)
# print(companies)

# problem 11
# data_dict = {
#     'company': ['Amazon', 'Microsoft', 'Facebook'],
#     'price': [2375.00, 178.6, 179.2],
#     'ticker': ['AMZN.US', 'MSFT.US', 'FB.US']
# }
#
# companies = pd.DataFrame(data=data_dict)
# companies = companies.set_index('company')
# print(companies)


# problem 12
# date_range = pd.date_range(start='2023-01-01', periods=31)
# print(date_range)

# problem 13


""" deleting >id from gff.3 files"""
# file=open('21-1.gff3.txt')
# a=open('my21.txt','a')
# for i in file:
#     if i.startswith(">"):
#         del(i)
#     else:
#         a.write(i)

"""sorting files have numbers listing inside the file using files"""
# file=open('newsort.txt')
# z=file.readlines()
# # x=z.sort()
# sorted(z)
# with open('sort.txt','w') as f:
#     for i in sorted(z):
#         f.write(i)


"""sorting files have numbers listing inside the file using pandas and write it inside excel file"""
# file=open('newsort.txt')
# z=file.readlines()
# r=[]
# for i in z:
#     r.append(i.split(' '))
#
# x=pd.DataFrame(data=r,columns=['num'+str(i) for i in range(20)])
# y=x.sort_values(by=['num0'])
# print(y)
# y.to_excel('file1.xlsx')


# def fibonacci():
#     num = int(input("How many numbers that generates?:")) #taking user input here with input function
#     # the value of iterator here is 1
#     i = 1
#     if num == 0: #if num == 0 which means that if number will = to 0 then
#         fib = [] #this string will be printed
#     elif num == 1: #if number == 1 so that will start from 1 and so on
#         fib = [1]
#     elif num == 2:
#         fib = [1,1]
#     elif num > 2:
#         fib = [1,1]
#         while i < (num - 1):
#             fib.append(fib[i] + fib[i-1]) #fibonacci logic
#             i += 1
#     return fib
# print(fibonacci()) #printing fibonacci funtion here
# input()

