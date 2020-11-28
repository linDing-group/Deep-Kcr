"""
在使用的时候，同种特征的数值变换（如1-space和2-space等）是在类(class)的层面；而进行序列的输入的时候都是在main()中；每个
函数按照提取到的特征的顺序定义了特征的名字，都在名为get_feature_name的函数中;对于AAindex这种特征来说，其特征维数和窗
口大小有关，所以在class的层面要给出序列！！
"""


class K_space:  # 注意，在使用get_feature_name的时候是默认输出0-k的名称，而在计算的时候只计算对应的那种k!!!
    def __init__(self, k):
        self.k = k
        self.feature_name = []
        self.AA_list_sort = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    def get_feature_name(self):
        for element in range(self.k + 1):
            element = str(element)
            prefix = 'Kspace_%s_space_' % element
            for i1 in self.AA_list_sort:
                for j1 in self.AA_list_sort:
                    self.feature_name.append(prefix + i1 + j1)
        return self.feature_name

    def main(self, sequence):
        two_mer_name = []
        for i2 in self.AA_list_sort:
            for j2 in self.AA_list_sort:
                combination = i2 + j2
                two_mer_name.append(combination)
        two_mer_dict = {}  # 定义一个2mer组成的字典用于计数
        for item in two_mer_name:
            two_mer_dict[item] = 0
        for i, v in enumerate(sequence):
            if i + 2 + self.k > len(sequence):  # 这里的思想是：i的索引加上1等于实际的位置，再加上k+1就是与之配对的那个（后面个）AA位置，看其超出序列的总长度没有
                break
            else:
                now_k_mer = sequence[i:i + 2 + self.k:self.k + 1]  # 后边界是取不到的，所以要定义为i+k+1
                if now_k_mer in two_mer_dict:  # 这里加一个判断是因为在补X的时候会有问题出现
                    two_mer_dict[now_k_mer] += 1
        feature_vector = []  # 这里可以直接加标签，通常以[]空的放置，[0]表示负样本标签,[1]表示正样本标签
        for key, value in two_mer_dict.items():
            frequency = value / (len(sequence) - 1 - self.k)
            feature_vector.append(frequency)
        return feature_vector


class CTD:
    def __init__(self):
        pass

# 首先进行转化:Polar(P);Neutral(N);Hydrophobicity(H);
# 转化为：Hydrophobicity
    def _transform1(self, sequence1):
        char_1 = ''
        class1 = 'RKEDQN'
        class2 = 'GASTPHY'
        class3 = 'CLVIMFW'
        for element1 in sequence1:
            if element1 in class1:
                char_1 += 'P'
            elif element1 in class2:
                char_1 += 'N'
            elif element1 in class3:
                char_1 += 'H'
            else:
                char_1 += 'X'
        return char_1

# 转化为：Normalized van der Waals volume
    def _transform2(self, sequence1):
        char_2 = ''
        class1 = 'GASTPDC'
        class2 = 'NVEQIL'
        class3 = 'MHKFRYW'
        for element2 in sequence1:
            if element2 in class1:
                char_2 += 'P'
            elif element2 in class2:
                char_2 += 'N'
            elif element2 in class3:
                char_2 += 'H'
            else:
                char_2 += 'X'
        return char_2

# 转化为：Polarity
    def _transform3(self, sequence1):
        char_3 = ''
        class1 = 'LIFWCMVY'
        class2 = 'PATGS'
        class3 = 'HQRKNED'
        for element3 in sequence1:
            if element3 in class1:
                char_3 += 'P'
            elif element3 in class2:
                char_3 += 'N'
            elif element3 in class3:
                char_3 += 'H'
            else:
                char_3 += 'X'
        return char_3

# 转化为：Polarizability
    def _transform4(self, sequence1):
        char_4 = ''
        class1 = 'GASDT'
        class2 = 'CPNVEQIL'
        class3 = 'KMHFRYW'
        for element4 in sequence1:
            if element4 in class1:
                char_4 += 'P'
            elif element4 in class2:
                char_4 += 'N'
            elif element4 in class3:
                char_4 += 'H'
            else:
                char_4 += 'X'
        return char_4

# 转化为：Charge
    def _transform5(self, sequence1):
        char_5 = ''
        class1 = 'KR'
        class2 = 'ANCQGHILMFPSTWYV'
        class3 = 'DE'
        for element5 in sequence1:
            if element5 in class1:
                char_5 += 'P'
            elif element5 in class2:
                char_5 += 'N'
            elif element5 in class3:
                char_5 += 'H'
            else:
                char_5 += 'X'
        return char_5

# 转化为：Secondary structure
    def _transform6(self, sequence1):
        char_6 = ''
        class1 = 'EALMQKRH'
        class2 = 'VIYCWFT'
        class3 = 'GNPSD'
        for element6 in sequence1:
            if element6 in class1:
                char_6 += 'P'
            elif element6 in class2:
                char_6 += 'N'
            elif element6 in class3:
                char_6 += 'H'
            else:
                char_6 += 'X'
        return char_6

# 转化为：Solvent accessibility
    def _transform7(self, sequence1):
        char_7 = ''
        class1 = 'ALFCGIVW'
        class2 = 'PKQEND'
        class3 = 'MRSTHY'
        for element7 in sequence1:
            if element7 in class1:
                char_7 += 'P'
            elif element7 in class2:
                char_7 += 'N'
            elif element7 in class3:
                char_7 += 'H'
            else:
                char_7 += 'X'
        return char_7

    def _computing_CTD(self, PNH_sequence):
        CTD_feature_vector = []
        length_sequence = len(PNH_sequence)
        # 计算C
        number_P = PNH_sequence.count('P')
        number_N = PNH_sequence.count('N')
        number_H = PNH_sequence.count('H')
        C_P = number_P / length_sequence
        C_N = number_N / length_sequence
        C_H = number_H / length_sequence
        CTD_feature_vector.append(C_P)
        CTD_feature_vector.append(C_N)
        CTD_feature_vector.append(C_H)
        # 计算T
        number_PN = PNH_sequence.count('PN')
        number_NP = PNH_sequence.count('NP')
        number_PN_NP = number_PN + number_NP
        number_PH = PNH_sequence.count('PH')
        number_HP = PNH_sequence.count('HP')
        number_PH_HP = number_PH + number_HP
        number_NH = PNH_sequence.count('NH')
        number_HN = PNH_sequence.count('HN')
        number_NH_HN = number_NH + number_HN
        T_PN_NP = number_PN_NP / (length_sequence - 1)
        T_PH_HP = number_PH_HP / (length_sequence - 1)
        T_NH_HN = number_NH_HN / (length_sequence - 1)
        CTD_feature_vector.append(T_PN_NP)
        CTD_feature_vector.append(T_PH_HP)
        CTD_feature_vector.append(T_NH_HN)
        # 计算D
        P_list = []
        N_list = []
        H_list = []
        for i, v in enumerate(PNH_sequence):
            if v == 'P':
                P_list.append(i)
            elif v == 'N':
                N_list.append(i)
            elif v == 'H':
                H_list.append(i)
        if P_list:
            first_P = (P_list[0] + 1) / length_sequence
            per25_P = (P_list[int((len(P_list) * 0.25)) - 1] + 1) / length_sequence
            per50_P = (P_list[int((len(P_list) * 0.5)) - 1] + 1) / length_sequence
            per75_P = (P_list[int((len(P_list) * 0.75)) - 1] + 1) / length_sequence
            per100_P = (P_list[-1] + 1) / length_sequence
        else:
            first_P = 0
            per25_P = 0
            per50_P = 0
            per75_P = 0
            per100_P = 0
        if N_list:
            first_N = (N_list[0] + 1) / length_sequence
            per25_N = (N_list[int((len(N_list) * 0.25)) - 1] + 1) / length_sequence
            per50_N = (N_list[int((len(N_list) * 0.5)) - 1] + 1) / length_sequence
            per75_N = (N_list[int((len(N_list) * 0.75)) - 1] + 1) / length_sequence
            per100_N = (N_list[-1] + 1) / length_sequence
        else:
            first_N = 0
            per25_N = 0
            per50_N = 0
            per75_N = 0
            per100_N = 0
        if H_list:
            first_H = (H_list[0] + 1) / length_sequence
            per25_H = (H_list[int((len(H_list) * 0.25)) - 1] + 1) / length_sequence
            per50_H = (H_list[int((len(H_list) * 0.5)) - 1] + 1) / length_sequence
            per75_H = (H_list[int((len(H_list) * 0.75)) - 1] + 1) / length_sequence
            per100_H = (H_list[-1] + 1) / length_sequence
        else:
            first_H = 0
            per25_H = 0
            per50_H = 0
            per75_H = 0
            per100_H = 0
        CTD_feature_vector.append(first_P)
        CTD_feature_vector.append(per25_P)
        CTD_feature_vector.append(per50_P)
        CTD_feature_vector.append(per75_P)
        CTD_feature_vector.append(per100_P)
        CTD_feature_vector.append(first_N)
        CTD_feature_vector.append(per25_N)
        CTD_feature_vector.append(per50_N)
        CTD_feature_vector.append(per75_N)
        CTD_feature_vector.append(per100_N)
        CTD_feature_vector.append(first_H)
        CTD_feature_vector.append(per25_H)
        CTD_feature_vector.append(per50_H)
        CTD_feature_vector.append(per75_H)
        CTD_feature_vector.append(per100_H)
        return CTD_feature_vector

    def get_feature_name(self):
        PNH_list = ['P', 'N', 'H']
        CTD_list = ['C', 'T', 'D']
        T_list = ['PN_NP', 'PH_HP', 'NH_HN']
        D_list = ['0', '25', '50', '75', '100']
        first_row_name = []
        for i in range(1, 8):
            for j in CTD_list:
                if j == 'C':
                    for k1 in PNH_list:
                        char = ('CTD_%a_' + j + '_' + k1) % i
                        first_row_name.append(char)
                elif j == 'T':
                    for k2 in T_list:
                        char = ('CTD_%a_' + j + '_' + k2) % i
                        first_row_name.append(char)
                elif j == 'D':
                    for k3 in PNH_list:
                        for l in D_list:
                            char = ('CTD_%a_' + j + '_' + k3 + '_' + l) % i
                            first_row_name.append(char)
        return first_row_name

    def main(self, sequence):
        t1 = self._transform1(sequence)
        feature1 = self._computing_CTD(t1)
        t2 = self._transform2(sequence)
        feature2 = self._computing_CTD(t2)
        t3 = self._transform3(sequence)
        feature3 = self._computing_CTD(t3)
        t4 = self._transform4(sequence)
        feature4 = self._computing_CTD(t4)
        t5 = self._transform5(sequence)
        feature5 = self._computing_CTD(t5)
        t6 = self._transform6(sequence)
        feature6 = self._computing_CTD(t6)
        t7 = self._transform7(sequence)
        feature7 = self._computing_CTD(t7)
        feature_total = feature1 + feature2 + feature3 + feature4 + feature5 + feature6 + feature7
        return feature_total


class EBGW:
    def __init__(self, number_sub_sequence):
        self.number_sub_sequence = number_sub_sequence

    def get_feature_name(self):
        feature_name = []
        for i in range(1, self.number_sub_sequence + 1):
            prefix = 'EBGW_%s_feature_' % i
            total_number_sub_sequence = 3 * i
            for j in range(1, total_number_sub_sequence + 1):
                featureName = prefix + str(j)
                feature_name.append(featureName)
        return feature_name

    def main(self, sequence):
        # The hydrophobic group
        c1 = 'AFGILMPVW'
        # The polar group
        c2 = 'CNQSTY'
        # The positively charged group
        c3 = 'KHR'
        # The negatively charged group
        c4 = 'DE'
        # definition of H1
        h1_bi_sequence = ''
        for element in sequence:
            if element in c1 or element in c2:
                h1_bi_sequence += str(1)
            else:
                h1_bi_sequence += str(0)

        # definition of H2
        h2_bi_sequence = ''
        for element in sequence:
            if element in c1 or element in c3:
                h2_bi_sequence += str(1)
            else:
                h2_bi_sequence += str(0)

        # definition of H3
        h3_bi_sequence = ''
        for element in sequence:
            if element in c1 or element in c4:
                h3_bi_sequence += str(1)
            else:
                h3_bi_sequence += str(0)

        feature = []
        feature_H1_list = []
        feature_H2_list = []
        feature_H3_list = []
        for i in range(1, self.number_sub_sequence + 1):
            i_length = round(i * len(sequence) / self.number_sub_sequence)  # definition length of the i-th sub-sequence
            i_h1_sub_bi_sequence = h1_bi_sequence[0:i_length]
            sum_i_h1 = i_h1_sub_bi_sequence.count('1')
            feature_h1_sub_sequence = sum_i_h1 / i_length
            feature_H1_list.append(feature_h1_sub_sequence)

            i_h2_sub_bi_sequence = h2_bi_sequence[0:i_length]
            sum_i_h2 = i_h2_sub_bi_sequence.count('1')
            feature_h2_sub_sequence = sum_i_h2 / i_length
            feature_H2_list.append(feature_h2_sub_sequence)

            i_h3_sub_bi_sequence = h3_bi_sequence[0:i_length]
            sum_i_h3 = i_h3_sub_bi_sequence.count('1')
            feature_h3_sub_sequence = sum_i_h3 / i_length
            feature_H3_list.append(feature_h3_sub_sequence)

        for element1 in feature_H1_list:
            feature.append(element1)
        for element2 in feature_H2_list:
            feature.append(element2)
        for element3 in feature_H3_list:
            feature.append(element3)
        return feature


class PWAA:  # 该种特征要取上下游等长的情况才行
    def __init__(self):
        self.AA_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    def get_feature_name(self):
        feature_name_list = []
        for element in self.AA_list:
            featureName = 'PWAA_feature_%s' % element
            feature_name_list.append(featureName)
        return feature_name_list

    def main(self, sequence):
        length_up_down = (len(sequence) - 1) / 2
        feature = []
        for aa_char in self.AA_list:
            sum_inter = 0
            if aa_char not in sequence:
                feature.append(0)
            else:
                for sequence_index, sequence_char in enumerate(sequence):
                    if sequence_char == aa_char:
                        j = sequence_index - length_up_down  # 这里10到时要改成上下游的那个L
                        sum_inter = sum_inter + (j + abs(j) / length_up_down)
                c = (1 / (length_up_down * (length_up_down + 1))) * sum_inter
                feature.append(c)
        return feature


# 使用AAindex特征需要提供AAidx数据文件
class AAindex:
    def __init__(self, length_sequence):
        import pandas as pd
        self.length_sequence = length_sequence
        self.obj = pd.read_csv('AAindex_12.csv')  # AAidx数据文件与该文件放在一起
        self.pro_name_list = self.obj['AccNo'].tolist()

    def get_feature_name(self):
        feature_name_list = []
        for i in range(1, self.length_sequence + 1):
            for element in self.pro_name_list:
                featureName = 'AAindex_pos%s_%s' % (i, element)
                feature_name_list.append(featureName)
        return feature_name_list

    def main(self, sequence):
        AA_list_sort = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                        'Y']
        AAindex_dict = {}
        for ele in AA_list_sort:
            AAindex_dict[ele] = self.obj[ele].tolist()
        AAindex_dict['X'] = [0] * 12
        feature = []
        for item in sequence:
            feature.extend(AAindex_dict[item])
        return feature


def main():
    f_positive = open('pos_train.fasta', 'r')
    f_negative = open('neg_train.fasta', 'r')
    f_total = open('5feature_idx12.csv', 'w')

    # 写第一行，特征名
    name_all_feature = []
    name_Kspace = K_space(5).get_feature_name()
    name_all_feature.extend(name_Kspace)
    # print(len(name_Kspace))
    name_AAindex = AAindex(31).get_feature_name()
    name_all_feature.extend(name_AAindex)
    # print(len(name_AAindex))
    name_CTD = CTD().get_feature_name()
    name_all_feature.extend(name_CTD)
    # print(len(name_CTD))
    name_EBGW = EBGW(5).get_feature_name()
    name_all_feature.extend(name_EBGW)
    # print(len(name_EBGW))
    name_PWAA = PWAA().get_feature_name()
    name_all_feature.extend(name_PWAA)
    # print(len(name_all_feature))
    # print('---------')
    # print(len(name_PWAA))
    f_total.write('class,')
    for ix, name in enumerate(name_all_feature):
        f_total.write(name)
        if ix != (len(name_all_feature) - 1):
            f_total.write(',')
        else:
            f_total.write('\n')

    # 写正样本特征
    for ele_pos in f_positive:
        if ele_pos[0] != '>':
            f_total.write('1')  # 标签！！！
            f_total.write(',')
            feature_total = []
            pos_sequence = ele_pos.strip()
            feature_0_space = K_space(0).main(pos_sequence)
            feature_1_space = K_space(1).main(pos_sequence)
            feature_2_space = K_space(2).main(pos_sequence)
            feature_3_space = K_space(3).main(pos_sequence)
            feature_4_space = K_space(4).main(pos_sequence)
            feature_5_space = K_space(5).main(pos_sequence)
            feature_total.extend(feature_0_space)
            feature_total.extend(feature_1_space)
            feature_total.extend(feature_2_space)
            feature_total.extend(feature_3_space)
            feature_total.extend(feature_4_space)
            feature_total.extend(feature_5_space)
            feature_AAindex = AAindex(31).main(pos_sequence)
            feature_total.extend(feature_AAindex)
            feature_CTD = CTD().main(pos_sequence)
            feature_total.extend(feature_CTD)
            feature_1_EBGW = EBGW(1).main(pos_sequence)
            feature_2_EBGW = EBGW(2).main(pos_sequence)
            feature_3_EBGW = EBGW(3).main(pos_sequence)
            feature_4_EBGW = EBGW(4).main(pos_sequence)
            feature_5_EBGW = EBGW(5).main(pos_sequence)
            feature_total.extend(feature_1_EBGW)
            feature_total.extend(feature_2_EBGW)
            feature_total.extend(feature_3_EBGW)
            feature_total.extend(feature_4_EBGW)
            feature_total.extend(feature_5_EBGW)
            feature_PWAA = PWAA().main(pos_sequence)
            feature_total.extend(feature_PWAA)
            for ix1, pos_feature in enumerate(feature_total):
                f_total.write(str(pos_feature))
                if ix1 != (len(feature_total) - 1):
                    f_total.write(',')
                else:
                    f_total.write('\n')

    # 写负样本特征
    for ele_neg in f_negative:
        if ele_neg[0] != '>':
            f_total.write('0')  # 标签！！！
            f_total.write(',')
            feature_total = []
            neg_sequence = ele_neg.strip()
            feature_0_space = K_space(0).main(neg_sequence)
            feature_1_space = K_space(1).main(neg_sequence)
            feature_2_space = K_space(2).main(neg_sequence)
            feature_3_space = K_space(3).main(neg_sequence)
            feature_4_space = K_space(4).main(neg_sequence)
            feature_5_space = K_space(5).main(neg_sequence)
            feature_total.extend(feature_0_space)
            feature_total.extend(feature_1_space)
            feature_total.extend(feature_2_space)
            feature_total.extend(feature_3_space)
            feature_total.extend(feature_4_space)
            feature_total.extend(feature_5_space)
            feature_AAindex = AAindex(31).main(neg_sequence)
            feature_total.extend(feature_AAindex)
            feature_CTD = CTD().main(neg_sequence)
            feature_total.extend(feature_CTD)
            feature_1_EBGW = EBGW(1).main(neg_sequence)
            feature_2_EBGW = EBGW(2).main(neg_sequence)
            feature_3_EBGW = EBGW(3).main(neg_sequence)
            feature_4_EBGW = EBGW(4).main(neg_sequence)
            feature_5_EBGW = EBGW(5).main(neg_sequence)
            feature_total.extend(feature_1_EBGW)
            feature_total.extend(feature_2_EBGW)
            feature_total.extend(feature_3_EBGW)
            feature_total.extend(feature_4_EBGW)
            feature_total.extend(feature_5_EBGW)
            feature_PWAA = PWAA().main(neg_sequence)
            feature_total.extend(feature_PWAA)
            for ix2, neg_feature in enumerate(feature_total):
                f_total.write(str(neg_feature))
                if ix2 != (len(feature_total) - 1):
                    f_total.write(',')
                else:
                    f_total.write('\n')


if __name__ == '__main__':
    main()
