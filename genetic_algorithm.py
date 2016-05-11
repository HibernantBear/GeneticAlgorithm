# coding=utf-8
import random

__author__ = "HibernantBear"


class GeneticAlgorithm(object):
    """
    遗传算法
    """

    def __init__(self):
        self.individuals = []  # 个体的集合
        self.individual_number = 0
        self.best_individual_ever = None  # 史上最好的个体
        self.cross_func = self._default_cross_func  # 交叉基因的函数，默认为_default_cross_func
        self.select_func = self._default_select_func  # 选择可以繁育下一代个体的函数，默认为_default_select_func
        self.fitting_func = None  # 适度计算函数，遗传算法的核心

    @staticmethod
    def _default_cross_func(parent1, parent2):
        """
        默认的基因交叉算法
        随机选取一个基因位置折断交换两个父本的基因
        :param parent1: 父本1
        :param parent2: 父本2
        :return: [子代1，子代2]
        """
        genes_num = len(parent1.DNA)
        break_point = random.randint(1, genes_num - 1)  # 基因折断点
        # 生成新个体并且赋予基因
        individual1 = Individual()
        individual1.dna_encode(*parent1.DNA[:-genes_num + break_point], *parent2.DNA[break_point:])
        individual1.set_variation_rate(parent1.variation_rate)
        individual2 = Individual()
        individual2.dna_encode(*parent2.DNA[:-genes_num + break_point], *parent1.DNA[break_point:])
        individual2.set_variation_rate(parent1.variation_rate)
        return [individual1, individual2]

    @staticmethod
    def _default_select_func(individual_all):
        """
        默认的选择函数
        舍弃群体里面所有适度值低的后半部分个体
        前半部分复制一份出来拥有两倍繁衍机会
        :param individual_all: 所有的个体
        :return: 选择后的群体
        """
        parents_num = len(individual_all)
        individual_all = individual_all[:-parents_num // 2]
        individual_all += copy.deepcopy(individual_all)
        return individual_all

    def individual_add(self, individual):
        self.individuals.append(individual)

    def get_best_individual_ever(self):
        return self.best_individual_ever

    def get_best_individual(self):
        return self.individuals[0]

    def get_worst_individual(self):
        return self.individuals[-1]

    def individual_variation(self, left_strong_one=True):
        for individual in self.individuals:
            if left_strong_one:
                individual_backup = copy.deepcopy(individual)
                individual.dna_variation()
                if self.fitting_func(*individual_backup.dna_decode()) > self.fitting_func(*individual.dna_decode()):
                    self.individuals[self.individuals.index(individual)] = individual_backup
            else:
                individual.dna_variation()

    def sort_by_fitting(self):
        fitting_result = {}  # 适度计算的结果
        assert self.fitting_func is not None  # 适度计算的函数存在
        for individual in self.individuals:
            fitting_result[individual] = self.fitting_func(*individual.dna_decode())  # 根据解码的DNA计算适度
        # 排序适度值
        self.individuals = [i[0] for i in sorted(fitting_result.items(), key=lambda x: x[1], reverse=True)]
        return max(fitting_result.items(), key=lambda x: x[1]), min(fitting_result.items(), key=lambda x: x[1])

    def update_best_individual_ever(self):
        self.sort_by_fitting()
        if self.fitting_func(*self.individuals[0].dna_decode()) > self.fitting_func(
                *self.best_individual_ever.dna_decode()):
            self.best_individual_ever = self.individuals[0]

    def limit_individual_number(self):
        while len(self.individuals) > self.individual_number:
            del self.individuals[-1]

    def _generate_offspring(self):
        rest_parents = self.select_func(self.individuals)  # 筛选
        # 打乱父本顺序为繁衍做准备
        random.shuffle(rest_parents)
        self.individuals = []  # 清空现在的群体记录
        for i in range(0, len(rest_parents), 2):
            self.individuals += self.cross_func(rest_parents[i], rest_parents[i + 1])
        self.individual_variation()  # 变异所有个体
        self.update_best_individual_ever()  # 更新最好的个体
        self.limit_individual_number()  # 限制群体数量

    def init_first_generation(self):
        self.sort_by_fitting()
        self.individual_number = len(self.individuals)
        self.best_individual_ever = self.individuals[0]

    def generate_offspring(self, n=0):
        self.init_first_generation()
        while True:
            self._generate_offspring()
            n -= 1
            if n == 0:
                break

    def offspring_generator(self, n=0, yield_all_individuals=False):
        self.init_first_generation()
        while True:
            self._generate_offspring()
            if yield_all_individuals:
                yield self.individuals
            else:
                yield self.get_best_individual()
            n -= 1
            if n == 0:
                break


class Individual(object):
    @staticmethod
    def _default_encode_func(*args):
        return args

    @staticmethod
    def _default_decode_func(dna):
        return dna

    @staticmethod
    def _default_variation_func(self, possible_genes):
        dna_length = len(self.DNA)
        variation_point = random.randint(1, dna_length)
        self.dna_encode(self.DNA[:-dna_length + variation_point - 1],
                        possible_genes[random.randint(0, len(possible_genes) - 1)],
                        self.DNA[variation_point:])

    encode_func = _default_encode_func
    decode_func = _default_decode_func
    variation_func = _default_variation_func
    possible_genes = []

    def __init__(self):
        self.DNA = None
        self.variation_rate = 0

    def dna_encode(self, *args):
        self.DNA = Individual.encode_func(*args)

    def dna_decode(self):
        return Individual.decode_func(self.DNA)

    def dna_variation(self):
        if random.random() < self.variation_rate:
            Individual.variation_func(self, Individual.possible_genes)

    def set_variation_rate(self, variation_rate):
        self.variation_rate = variation_rate


if __name__ == "__main__":
    import copy

    l1 = [1, 2, 3, 4, 5, 6, 7]
    l2 = [2, 5, 8, 9, 4, 6, 1]

    ga = GeneticAlgorithm()


    def fitting_func_custom(*args):
        tmp = 0
        for arg in args:
            tmp += arg
        return tmp


    ga.fitting_func = fitting_func_custom

    for i_test in range(4):
        individual_test = Individual()
        individual_test.dna_encode(l1[random.randint(0, 6)], l2[random.randint(0, 6)])
        individual_test.set_variation_rate(.5)
        ga.individual_add(individual_test)
    Individual.possible_genes = [l1] + [l2]


    def variation_func_custom(self, possible_genes):
        dna_length = len(self.DNA)
        variation_point = random.randint(1, dna_length)
        self.dna_encode(*self.DNA[:-dna_length + variation_point - 1],
                        possible_genes[variation_point - 1][random.randint(0, len(possible_genes) - 1)],
                        *self.DNA[variation_point:])


    Individual.variation_func = variation_func_custom

    for individuals_test in ga.offspring_generator(n=10, yield_all_individuals=True):
        s = ""
        for individual_test in individuals_test:
            s += str(individual_test.DNA)
        print(s + "\n")
