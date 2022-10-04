#    _____       __              
#   / ___/____  / /   _____  _____
#   \__ \/ __ \/ / | / / _ \/ ___/
#  ___/ / /_/ / /| |/ /  __/ /    
# /____/\____/_/ |___/\___/_/     
                                
from colorama import Fore, Back, Style, init
from prettytable import PrettyTable
import sympy as sp

class Solver(): 
    '''
    Класс реализации 
    
    Аргументы: 
        - погрешность -> e 
        - левый конец интервала -> a
        - правый конец интервала -> b 
    '''
    def __init__(self, e, a, b):
        init() # colorama
        self.e = e
        self.a = a
        self.b = b 
        self.x = sp.Symbol('x')
        self.func = 5 * self.x - 8 * sp.log(self.x) - 8
        #x^3 + 5*x^2+6*x
        self.func = sp.Pow(self.x,3)+5*sp.Pow(self.x,2)+6*self.x
        self.set_table_combine_method()
        self.set_table_simple_iteration_method()
        self.get_solv_combine_method()
        # self.get_solv_simple_iteration_method()
    
    def set_table_combine_method(self):
        '''
        Создание объекта класса PrettyTable. Установление полей для таблицы 
        
        `x0` - начальная точка для касательной
         
        `c` - начальная точка для хорды
        
        `x1` - точка пересечения касательной и оси x 
        
        `z` - точка пересечения хорды и оси x
        
        `r` - итерационная разность
        '''
        self.table_combine_method = PrettyTable() 
        self.table_combine_method.field_names = ['iteration','x0', 'c', 'x1','z','r']  

    
    def set_table_simple_iteration_method(self):
        '''
        Создание объекта класса PrettyTable. Установление полей для таблицы 
        '''
        self.table_simple_iteration_method = PrettyTable()
        self.table_simple_iteration_method.field_names = ['iteration','x0', 'c', 'x1','z','r'] # R - инетарционная разность 
    

    def get_start_data(self):
        '''
        Вывод входных данных
        '''
        print(f"{Fore.RED}Уравнение{Style.RESET_ALL}:    5 * x - 8 * ln(x) - 8 = 0")
        print(f"{Fore.RED}Погрешность{Style.RESET_ALL}:  {self.e} ")
        print(f"{Fore.RED}Интервал{Style.RESET_ALL}:     [{self.a} ; {self.b}] {Style.RESET_ALL}\n")     
    

    def get_about_method(self, flag: int):
        print(f"\n{Fore.GREEN}----------------------------------------{Style.RESET_ALL}")
        print(flag := 'Комбинированный метод' if (flag == 1) else 'Метод простой итерации') 
        print(f"{Fore.GREEN}----------------------------------------{Style.RESET_ALL}")
    

    def get_solv_combine_method(self):
        '''
        Суть метода - сужение интервала изоляции с двух сторон 
        
        На одной стороне сужение с помощью хорды на другой с и спользованием касательной          
        '''
        
        self.get_about_method(flag=1)
        self.get_start_data()
        # 5 * self.x - 8 * sp.log(self.x) - 8
        derivative_f = self.func.diff(self.x)
        derivative_f_2_order = derivative_f.diff(self.x)
        f_a = sp.Pow(self.a,3)+5*sp.Pow(self.a,2)+6*(self.a) # f(-2.5)
        f_b = sp.Pow(self.b,3)+5*sp.Pow(self.b,2)+6*(self.b) # f(-1.7)  x^3 + 5*x^2+6*x
        f_a_derivative = derivative_f.subs(self.x, self.a)
        f_b_derivative = derivative_f.subs(self.x, self.b)
        f_a_derivative_2_order = derivative_f_2_order.subs(self.x, self.a)
        f_b_derivative_2_order = derivative_f_2_order.subs(self.x, self.b)
        
        print(f"f`(x) = {derivative_f}\n"
              f"f``(x) = {derivative_f_2_order}\n"
              f"f(a) = {f_a}\n"
              f"f(b) = {f_b}\n"
              f"f'(a) = {f_a_derivative}\n"
              f"f'(b) = {f_b_derivative}\n"
              f"f''(a) = {f_a_derivative_2_order}\n"
              f"f''(b) = {f_b_derivative_2_order}\n")

        if (f_a * f_a_derivative_2_order > 0): 
            print((f'Неподвижна в точке a =>  {Fore.GREEN}по недостатку методом касательной \n\t\t\t по избытку методом хорд{Style.RESET_ALL}'))
        if (f_b * f_a_derivative_2_order > 0): 
            print(f'Неподвижна в точке b =>  {Fore.GREEN}по недостатку методом хорд \n\t\t\t по избытку методом касательных{Style.RESET_ALL}')        
        
        iteration = 1
        x0 = self.b
        c = self.a 
        r = 1 # заглушка 
        while (self.e <= r):
            x1 = x0-(sp.Pow(x0,3) + 5 * sp.Pow(x0, 2) + 6 * x0) / (3 * sp.Pow(x0, 2) + 10 * x0 + 6) # проводим касательную
            z = c-(sp.Pow(c,3) + 5 * sp.Pow(c, 2) + 6 * c) * (x0 - c) / ((sp.Pow(x0, 3) + 5 * sp.Pow(x0, 2) + 6 * x0) - (sp.Pow(c, 3) + 5 * sp.Pow(c, 2) + 6 * c)) # проводим хорду 
            r = abs(x1 - z) # итерационная разность 
            x0 = x1
            c = z              
            self.table_combine_method.add_row([iteration, x0, c, x1, z, r])
            iteration += 1 
      
        print(self.table_combine_method)
        print(f'Ответ получен в {iteration - 1} итерации\nОтвет: {Fore.GREEN}{x1}{Style.RESET_ALL}')

    def get_solv_simple_iteration_method(self): 
        self.get_about_method(flag=0)
        self.get_start_data()


if __name__ == "__main__": 
    solv = Solver(1e-5, -2.5, -1.7)
        