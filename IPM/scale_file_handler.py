### Scale file handler v1
### Jpevey 4/24/20
import collections
import os


class scale_file_handler:

    def __init__(self):
        ###### Fix these imports

        print("Let's play with some scale files!")
        self.data_dict = collections.OrderedDict()
        self.data_list = []

    ### Functions Related to materials

    ### Materials from Interior_Point_Materials_()_Cd.xls
    def build_default_material_dicts(self):

        materials_original = collections.OrderedDict()
        materials_original['fuel'] = collections.OrderedDict()
        materials_original['moderator'] = collections.OrderedDict()
        materials_original['poison'] = collections.OrderedDict()
        materials_original['stainless_steel_304'] = collections.OrderedDict()
        materials_original['void'] = collections.OrderedDict()

        materials_original['fuel']['u-234'] = 5.99675E-06
        materials_original['fuel']['u-235'] = 6.73792E-04
        materials_original['fuel']['u-236'] = 3.09945E-06
        materials_original['fuel']['u-238'] = 2.72872E-03
        materials_original['fuel']['n-14'] = 3.41161E-03
        materials_original['fuel']['si-28'] = 1.18202E-02
        materials_original['fuel']['si-29'] = 6.00478E-04
        materials_original['fuel']['si-30'] = 3.96303E-04
        materials_original['fuel']['c'] = 2.36746E-02
        materials_original['fuel']['he-3'] = 1.83606E-10
        materials_original['fuel']['he-4'] = 1.83606E-04

        materials_original['moderator']['y-89'] = 1.41215E-02
        materials_original['moderator']['h-1'] = 2.40037E-02
        materials_original['moderator']['h-2'] = 2.76075E-06
        materials_original['moderator']['si-28'] = 1.75075E-02
        materials_original['moderator']['si-29'] = 8.89400E-04
        materials_original['moderator']['si-30'] = 5.86983E-04
        materials_original['moderator']['c'] = 1.89838E-02

        materials_original['poison']['cd-106'] = 5.79177E-04
        materials_original['poison']['cd-108'] = 4.12372E-04
        materials_original['poison']['cd-110'] = 5.77787E-03
        materials_original['poison']['cd-111'] = 5.93074E-03
        materials_original['poison']['cd-112'] = 1.11711E-02
        materials_original['poison']['cd-113'] = 5.66665E-03
        materials_original['poison']['cd-114'] = 1.33210E-02
        materials_original['poison']['cd-116'] = 3.47968E-03

        materials_original['stainless_steel_304']['c'] = 6.02E-05
        materials_original['stainless_steel_304']['si-28'] = 7.94E-04
        materials_original['stainless_steel_304']['si-29'] = 3.91E-05
        materials_original['stainless_steel_304']['si-30'] = 2.49E-05
        materials_original['stainless_steel_304']['p-31'] = 3.58E-05
        materials_original['stainless_steel_304']['s-32'] = 2.15E-05
        materials_original['stainless_steel_304']['s-33'] = 1.64E-07
        materials_original['stainless_steel_304']['s-34'] = 9.04E-07
        materials_original['stainless_steel_304']['s-36'] = 2.01E-09
        materials_original['stainless_steel_304']['cr-50'] = 7.96E-04
        materials_original['stainless_steel_304']['cr-52'] = 1.48E-02
        materials_original['stainless_steel_304']['cr-53'] = 1.64E-03
        materials_original['stainless_steel_304']['cr-54'] = 4.01E-04
        materials_original['stainless_steel_304']['mn-55'] = 8.77E-04
        materials_original['stainless_steel_304']['fe-54'] = 3.63E-03
        materials_original['stainless_steel_304']['fe-56'] = 5.49E-02
        materials_original['stainless_steel_304']['fe-57'] = 1.25E-03
        materials_original['stainless_steel_304']['fe-58'] = 1.62E-04
        materials_original['stainless_steel_304']['ni-58'] = 5.66E-03
        materials_original['stainless_steel_304']['ni-60'] = 2.11E-03
        materials_original['stainless_steel_304']['ni-61'] = 9.01E-05
        materials_original['stainless_steel_304']['ni-62'] = 2.83E-04
        materials_original['stainless_steel_304']['ni-64'] = 6.98E-05

        materials_original['void']["""\'"""] = 0
        return materials_original


    def build_material_dictionaries(self, default_materials_definition_from_options):
        list_of_material_dictionaries = []
        materials = self.build_default_material_dicts()
        for material_definition in default_materials_definition_from_options:
            ### If it is a mixture of two or more materials:
            if ":" in material_definition:
                mat_dif_split = material_definition.split(":")
                local_material_list = mat_dif_split[0].split("/")
                local_material_ratio = mat_dif_split[1].split("/")

                if len(local_material_list) > 2:
                    print("Unable to combine more than 2 materials into a default material... fix this")
                    exit()

                mat_dict = self.combine_material_dicts(materials[local_material_list[0]],
                                                       materials[local_material_list[1]],
                                                       int(local_material_ratio[0]))
            else:
                mat_dict = materials[material_definition]

            list_of_material_dictionaries.append(mat_dict)
        return list_of_material_dictionaries

    def perturb_dict(self, data_dict, perturbation=0):
        perturbation = perturbation / 100

        new_dict = collections.OrderedDict()
        for key in data_dict:
            new_dict[key] = data_dict[key] * perturbation

        return new_dict

    def build_scale_material_string(self, material_dictionary, material_number, temperature):
        isotope_string = ""
        for isotope in material_dictionary:
            isotope_nuclear_density = str(material_dictionary[isotope])
            temperature = str(temperature)

            material_number = str(material_number)

            isotope_string += str(
                isotope) + " " + material_number + " 0 " + isotope_nuclear_density + " " + temperature + " end\n"

        return isotope_string

    def combine_material_dicts(self, material_1_dict, material_2_dict, beta):
        new_material_dict = collections.OrderedDict()
        # print(material_1_dict)
        if beta >= 1:
            beta = beta / 100

        material_1_beta = beta
        material_2_beta = (1 - beta)

        ### Creating list of all materials
        list_of_materials = []
        for material in material_1_dict:
            if material not in list_of_materials:
                list_of_materials.append(material)
                # print(material)
        for material in material_2_dict:
            if material not in list_of_materials:
                list_of_materials.append(material)
                # print(material)
        # print(list_of_materials)
        for isotope in list_of_materials:
            mat_1_val = 0
            mat_2_val = 0

            try:
                mat_1_val = material_1_dict[isotope] * material_1_beta
            except:
                pass

            try:
                mat_2_val = material_2_dict[isotope] * material_2_beta
            except:
                pass

            new_material_dict[isotope] = mat_1_val + mat_2_val
            # print(isotope, new_material_dict[isotope])

        return new_material_dict

    ### Dealing with getting data out of files
    def get_keff_for_all_scale_outputs(self, write_out_data_dict_bool=True, all_output_files=True):
        if all_output_files:
            for file in os.listdir("."):
                if file.endswith('.out'):
                    print("Found output file:", file)
                    self.get_keff_and_uncertainty(file)
        if write_out_data_dict_bool:
            self.writeout_data_dict(output_file_string)

    def get_keff_and_uncertainty(self, file_name):
        if 'keff' not in self.data_list:
            self.data_list.append('keff')
        if 'keff_uncertainty' not in self.data_list:
            self.data_list.append('keff_uncertainty')

        if file_name not in self.data_dict:
            self.data_dict[file_name] = collections.OrderedDict()

        output_file = open(file_name, 'r')
        for line in output_file:
            if "best estimate system k-eff" in line:
                line_split = line.split("                               ")
                line_split_2 = line_split[1].split(" + or -")

                self.data_dict[file_name]['keff'] = line_split_2[0].strip()

                line_split_3 = line_split_2[1].strip().split(" ")
                self.data_dict[file_name]['keff_uncertainty'] = line_split_3[0].strip()

        print("found keff:", self.data_dict[file_name]['keff'], " +/- ", self.data_dict[file_name]['keff_uncertainty'])
        return self.data_dict[file_name]['keff'], self.data_dict[file_name]['keff_uncertainty']

    def writeout_data_dict(self, output_file_string):
        print("Writing out file: ", output_file_string)
        output_file = open(output_file_string, "w")

        header = ""
        for file in self.data_dict:
            header += "," + file
        output_file.write(header + "\n")

        for data_type in self.data_list:
            write_string = str(data_type)
            for file in self.data_dict:
                write_string += "," + self.data_dict[file][data_type]
            output_file.write(write_string + "\n")
        output_file.close()

    ### Tsunami File function section
    def parse_sdf_file_into_dict(self, tsunami_file_string):
        tsunami_file = open(tsunami_file_string, 'r')
        in_data = False
        data_dict = collections.OrderedDict()
        for line in tsunami_file:
            line = line.strip()
            if line == "0.000000E+00  0.000000E+00      0      0":
                continue
            if 'total' in line:
                in_data = True
                line = line.strip()

                line_split = line.split('total')
                isotope = line_split[0].strip()

                in_data_count = 0
                continue

            if in_data:
                # print(line.strip())
                if in_data_count == 0:
                    line_split = line.split(' ')
                    material = line_split[0]

                if in_data_count == 1:
                    line_split = line.split('  ')
                    sensitivity = line_split[0]
                    uncert = line_split[1]

                in_data_count += 1
                if in_data_count == 2:
                    in_data = False
                    # print(material, isotope, sensitivity, uncert)
                    if material not in data_dict:
                        data_dict[material] = collections.OrderedDict()
                    data_dict[material][isotope] = collections.OrderedDict()
                    data_dict[material][isotope]['sensitivity'] = sensitivity
                    data_dict[material][isotope]['uncertainty'] = uncert
        return data_dict

    def build_scale_input_from_beta(self,
                                    material_betas,
                                    material_1,
                                    material_2,
                                    template_file_string,
                                    flag,
                                    flag_replacement_string='replace',
                                    temperature=300,
                                    material_count_offset=1,
                                    file_name_flag='default_'):
        material_list = []
        for beta in material_betas:
            material_list.append(self.combine_material_dicts(material_1, material_2, beta))

        material_string_list = []
        for count, material in enumerate(material_list):
            material_string_list.append(
                self.build_scale_material_string(material, count + material_count_offset, temperature))

        ### Making list of keys
        flag_list = []
        for x in range(len(material_string_list)):
            flag_list.append(flag.replace(flag_replacement_string, str(x)))

        material_dict = self.make_data_dict(flag_list, material_string_list)

        self.create_scale_input_given_target_dict(template_file_string, file_name_flag, material_dict)

    def writeout_total_sensitivity_per_material_per_isotope(self, data_dict, output_file_string):
        output_file = open(output_file_string, 'w')
        for material in data_dict:
            for isotope in data_dict[material]:
                output_file.write(material + "," + isotope + "," + data_dict[material][isotope]['sensitivity'] + "," +
                                  data_dict[material][isotope]['uncertainty'] + "\n")

    def writeout_total_sensitivity_per_isotope_per_material_one_table(self, data_dict, output_file_string, tables):
        output_file = open(output_file_string, 'w')
        tables = ['sensitivity', 'uncertainty']

        ### building headers
        example_isotope = ""
        example_material = ""
        material_header = ""
        for material in data_dict:
            if example_material == "":
                example_material = material
            material_header += "," + material
            # print(material)
            for isotope in data_dict[material]:
                if example_isotope == "":
                    example_isotope = isotope

        # print(example_material, example_isotope)

        for table_val in tables:
            output_file.write(table_val + "\n")
            output_file.write(material_header + "\n")
            for isotope in data_dict[example_material]:

                write_string = ""
                for material in data_dict:
                    write_string += data_dict[material][isotope][table_val] + ","
                output_file.write(isotope + "," + write_string + "\n")
            output_file.write("\n")

        output_file.close()

    ### Given a list of values, return a new list with those multiplied
    ### by a pertubation. This perturbation is divided by 100 (i.e. +5% would be 5)
    def return_perturbed_materials(self, list_of_materials, perturbation):
        perturbation = float(perturbation) / 100
        new_list = []
        for item in list_of_materials:
            new_list.append(float(item) + float(item) * perturbation)
        return new_list

    ### This function turns a list of keys and a list of data into an
    ### ordereddict()
    def make_data_dict(self, keys_list, data_list):
        data_dict = collections.OrderedDict()
        for item, data_value in zip(keys_list, data_list):
            data_dict[item] = data_value
        return data_dict

    ### Takes a flag, a keyword, and a number of flags to return and returns a list of the flags
    def return_target_list(self, flag, flag_keyword, number_of_targets):
        target_list = []
        for _ in range(number_of_targets):
            target_list.append(flag.replace(flag_keyword, str(_)))
        return target_list

    ### This function takes a dictionary of keys/data and opens a template file,
    ### replacing the keys found with the data. Returns the name of the file created.
    def create_scale_input_given_target_dict(self, template_file_string, flag_, target_dict):
        run_file_list = []
        template_file = open(template_file_string, 'r')
        output_filename = flag_ + ".inp"
        output_file = open(output_filename, 'w')

        for line in template_file:
            for target_string in target_dict:
                if target_string in line:
                    line = line.replace(str(target_string), str(target_dict[target_string]))
            output_file.write(line)
        output_file.close()
        template_file.close()

        return output_filename
