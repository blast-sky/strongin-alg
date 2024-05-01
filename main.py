# noinspection PyUnresolvedReferences
from math import *
from tkinter import *

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

STRONGIN_CONSTANT = "Стронгина"
PIYAVSKIY_CONSTANT = "Пиявского"

SMALL_PADDING = {"padx": 5, "pady": 5}
PADDING = {"padx": 10, "pady": 10}
INPUT_PADDING = {"ipadx": 2, "ipady": 2}

IDENTITY_PARAMETER_COMPUTER_METHOD = "IDENTITY_METHOD"
STRONGIN_PARAMETER_COMPUTER_METHOD = "STRONGIN_METHOD"


class StronginMethodApp:
    def __init__(self):
        self.root = Tk()
        self.root.title("Поиск глобального минимума функции")
        self.root.geometry("1200x600")
        self.root.resizable(False, False)

        self.plot_widget: Canvas | None = None

        self.configuration_frame = Frame(self.root)

        function_frame = self.create_function_frame()
        function_frame.grid(row=0, column=0, columnspan=2, **PADDING)

        border_frame = self.create_border_frame()
        border_frame.grid(row=1, column=0, columnspan=2, **PADDING)

        parameter_frame = self.create_parameter_frame()
        parameter_frame.grid(row=2, column=0, columnspan=2, padx=10, pady=20)

        parameter_method_label = Label(self.configuration_frame, text="Подсчет характеристики:")
        parameter_method_label.grid(row=5, column=0, columnspan=2)

        self.selected_parameter_computer_method = StringVar(value=IDENTITY_PARAMETER_COMPUTER_METHOD)

        parameter_computer_methods = [
            ("Фиксированное значение", IDENTITY_PARAMETER_COMPUTER_METHOD),
            ("Адаптивное значение", STRONGIN_PARAMETER_COMPUTER_METHOD)
        ]

        current_row = 6
        for text, value in parameter_computer_methods:
            method_radio_button = Radiobutton(
                self.configuration_frame,
                text=text,
                variable=self.selected_parameter_computer_method,
                value=value
            )
            method_radio_button.grid(row=current_row, column=0, sticky=W, columnspan=2, padx=100)
            current_row += 1

        self.confirm_button = Button(self.configuration_frame, text="Начать", command=self.create_prot)
        self.confirm_button.grid(row=8, column=0, columnspan=2, pady=30, **INPUT_PADDING)

        result_frame = self.create_result_frame()
        result_frame.grid(row=9, column=0, columnspan=2)

        self.configuration_frame.pack(side=LEFT, fill=Y)

        self.root.mainloop()

    def create_parameter_frame(self):
        parameter_frame = Frame(self.configuration_frame)

        r_label = Label(parameter_frame, text="Параметр (r):")
        r_label.grid(row=2, column=0, **SMALL_PADDING, sticky=E)

        self.r_entry = Entry(parameter_frame)
        self.r_entry.insert(0, "2")
        self.r_entry.grid(row=2, column=1, **SMALL_PADDING)

        error_label = Label(parameter_frame, text="Погрешность:")
        error_label.grid(row=3, column=0, **SMALL_PADDING, sticky=E)

        self.error_entry = Entry(parameter_frame)
        self.error_entry.insert(0, "0.001")
        self.error_entry.grid(row=3, column=1)

        max_iteration_label = Label(parameter_frame, text="Максимальное число итераций:")
        max_iteration_label.grid(row=4, column=0, **SMALL_PADDING, sticky=E)

        self.max_iteration_entry = Entry(parameter_frame)
        self.max_iteration_entry.insert(0, "800")
        self.max_iteration_entry.grid(row=4, column=1)

        return parameter_frame

    def create_function_frame(self):
        function_frame = Frame(self.configuration_frame)

        function_label = Label(function_frame, text="F(x) = ")
        function_label.pack(side=LEFT)

        self.function_entry = Entry(function_frame, width=30)
        self.function_entry.insert(0, "2 * sin(3 * x) + 3 * cos(5 * x)")
        self.function_entry.pack(side=LEFT, **INPUT_PADDING)

        return function_frame

    def create_border_frame(self):
        border_frame = Frame(self.configuration_frame)

        border_label_1 = Label(border_frame, text="x ∈ [")
        border_label_1.pack(side=LEFT)

        self.left_border_entry = Entry(border_frame, width=5)
        self.left_border_entry.insert(0, "0")
        self.left_border_entry.pack(side=LEFT, **INPUT_PADDING)

        border_label_2 = Label(border_frame, text=" , ")
        border_label_2.pack(side=LEFT)

        self.right_border_entry = Entry(border_frame, width=5)
        self.right_border_entry.insert(0, "8")
        self.right_border_entry.pack(side=LEFT, **INPUT_PADDING)

        border_label_3 = Label(border_frame, text="]")
        border_label_3.pack(side=LEFT)

        return border_frame

    def create_result_frame(self):
        result_frame = Frame(self.configuration_frame)

        minimum_label = Label(result_frame, text="Найденный минимум:")
        minimum_label.grid(row=2, column=0, **SMALL_PADDING, sticky=E)

        self.minimum_entry = Entry(result_frame)
        self.minimum_entry.grid(row=2, column=1, **SMALL_PADDING)

        minimum_point_label = Label(result_frame, text="Точка минимума:")
        minimum_point_label.grid(row=3, column=0, **SMALL_PADDING, sticky=E)

        self.minimum_point_entry = Entry(result_frame)
        self.minimum_point_entry.grid(row=3, column=1)

        iteration_count_label = Label(result_frame, text="Число итераций:")
        iteration_count_label.grid(row=4, column=0, **SMALL_PADDING, sticky=E)

        self.iteration_count_entry = Entry(result_frame)
        self.iteration_count_entry.grid(row=4, column=1)

        return result_frame

    def create_prot(self):
        plot_configurator = PlotConfigurator(
            self.root,
            self.function_entry.get(),
            self.left_border_entry.get(),
            self.right_border_entry.get(),
            self.r_entry.get(),
            self.selected_parameter_computer_method.get(),
            self.max_iteration_entry.get(),
            self.error_entry.get()
        )

        if self.plot_widget is not None:
            self.plot_widget.pack_forget()

        self.plot_widget = plot_configurator.get_widget()
        self.plot_widget.pack(side=RIGHT, fill=BOTH, expand=True, padx=10, pady=10)

        point, minimum, iteration_count = plot_configurator.get_results()

        self.minimum_entry.delete(0, END)
        self.minimum_entry.insert(0, minimum)
        self.iteration_count_entry.delete(0, END)
        self.iteration_count_entry.insert(0, str(iteration_count))
        self.minimum_point_entry.delete(0, END)
        self.minimum_point_entry.insert(0, str(point))


class PlotConfigurator:
    def __init__(self, master, function_string, a, b, r, parameter_compute_method, max_iteration_count, error):
        self.master = master
        self.function_string = function_string
        self.a = float(a)
        self.b = float(b)
        self.r = float(r)
        self.parameter_compute_method = parameter_compute_method
        self.max_iteration_count = int(max_iteration_count)
        self.error = float(error)

        if self.parameter_compute_method == STRONGIN_PARAMETER_COMPUTER_METHOD:
            self.iteration_function = self.strongin_iteration
        else:
            self.iteration_function = self.piyavskiy_iteration

        self.points = [(self.a, self.compute_function(self.a), 1), (self.b, self.compute_function(self.b), 2)]

        self.iteration_count = 2

        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        # self.canvas_widget = self.canvas.get_tk_widget()
        # self.canvas_widget.pack(side=LEFT, fill=BOTH, expand=True, padx=10, pady=10)

        self.compute_points()
        self.plot_function()
        self.canvas.draw()

    def compute_points(self):
        while self.iteration_count < self.max_iteration_count and self.min_diff_in_points() > self.error:
            self.iteration_count += 1
            self.iterate()

    def get_widget(self):
        return self.canvas.get_tk_widget()

    def get_results(self):
        min_func_value_point = min(self.points, key=lambda x: x[1])
        return min_func_value_point[0], min_func_value_point[1], self.iteration_count

    def min_diff_in_points(self):
        min_diff = self.points[1][0] - self.points[0][0]
        for i in range(3, len(self.points)):
            min_diff = min(min_diff, self.points[i][0] - self.points[i - 1][0])
        return min_diff

    def plot_function(self):
        x = np.linspace(self.a, self.b, 1000)
        y = [self.compute_function(cur_x) for cur_x in x]
        self.ax.plot(x, y, label='Исходная функция')

        points_x = [point[0] for point in self.points]
        points_y = [point[1] for point in self.points]
        self.ax.scatter(points_x, points_y, color='red', label='Точки испытаний')
        # self.ax.plot(points_x, points_y, color='red', linestyle='dashed')

        for point in self.points:
            x, y, iteration = point
            # self.ax.text(x, y, f'{iteration}', fontsize=10, ha='right', va='bottom')

        self.ax.legend()

    def iterate(self):
        x_new, R_max = self.iteration_function()

        self.points.append((x_new, self.compute_function(x_new), self.iteration_count))

        self.points.sort(key=lambda point: point[0])

        # print(f"Итерация #{self.iteration_count}: x = {x_new:.6f}, z = {self.compute_function(x_new):.6f}, max R = {R_max:.3f}\n")

        # self.fig.clear()
        # self.ax = self.fig.add_subplot(111)
        # self.plot_function()
        #
        # self.canvas.draw()

    def compute_function(self, x):
        x = float(x)
        return eval(self.function_string, globals(), {"x": x})

    def strongin_iteration(self):
        M = []
        for i in range(len(self.points) - 1):
            x_i = self.points[i][0]
            x_i_plus_1 = self.points[i + 1][0]

            z_i = self.points[i][1]
            z_i_plus_1 = self.points[i + 1][1]

            M.append(abs(z_i - z_i_plus_1) / (x_i_plus_1 - x_i))
        # print(f"Iteration #{self.iteration_count} - M list: {M}")

        max_M = max(M)
        # print(f"Iteration #{self.iteration_count} - max M: {max_M}")

        m = self.r * max_M if max_M > 0 else 1

        R_list = []
        for i in range(len(self.points) - 1):
            R_list.append(self.piyavskiy_characteristic(self.points[i], self.points[i + 1], m))
        # print(f"Iteration #{self.iteration_count} - R list: {R_list}")
        max_R_index = R_list.index(max(R_list))

        interval_left_index = max_R_index
        interval_right_index = max_R_index + 1

        x_i = self.points[interval_left_index][0]
        x_i_plus_1 = self.points[interval_right_index][0]

        z_i = self.points[interval_left_index][1]
        z_i_plus_1 = self.points[interval_right_index][1]

        return 0.5 * (x_i + x_i_plus_1) - (z_i_plus_1 - z_i) / (2 * m), R_list[max_R_index]

    # def strongin_characteristic(self, left, right, m):
    #     return m * (right[0] - left[0]) + \
    #         (left[1] - right[1]) ** 2 / (m * (right[0] - left[0])) - \
    #         2 * (right[1] + left[1])

    def piyavskiy_iteration(self):
        R_list = []
        for i in range(len(self.points) - 1):
            R_list.append(self.piyavskiy_characteristic(self.points[i], self.points[i + 1], self.r))
        max_R_index = R_list.index(max(R_list))
        # print(f"Iteration #{self.iteration_count + 1} - R list: {R_list}")

        interval_left_index = max_R_index
        interval_right_index = max_R_index + 1

        x_i = self.points[interval_left_index][0]
        x_i_plus_1 = self.points[interval_right_index][0]

        z_i = self.points[interval_left_index][1]
        z_i_plus_1 = self.points[interval_right_index][1]

        return 0.5 * (x_i_plus_1 + x_i) - (z_i_plus_1 - z_i) / (2 * self.r), R_list[max_R_index]

    def piyavskiy_characteristic(self, left, right, r):
        return 0.5 * r * (right[0] - left[0]) - (right[1] + left[1]) / 2


if __name__ == "__main__":
    StronginMethodApp()
