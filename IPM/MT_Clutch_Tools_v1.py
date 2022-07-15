#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import socket
import time
import scale_file_handler
import collections


# In[ ]:


class MT_Clutch_Tools:
    def __init__(self, template_file = "", neutrons_per_generation = 0, skip_generations= 0, list_of_material_dictionaries = ""):
        self.neutrons_per_generation = neutrons_per_generation
        self.skip_generations = skip_generations
        self.template_file = template_file
        self.list_of_material_dictionaries = list_of_material_dictionaries

        self.sfh = scale_file_handler.scale_file_handler()

        print("Multithreaded Clutch tool")
        print("Scale settings:")
        print("neutrons per generation:", self.neutrons_per_generation)
        print("skip_generations:", self.skip_generations)
        print("clutch template file:", self.template_file)

        self.singlethreaded_scale_script = """#!/bin/bash

#PBS -q corei7
#PBS -V
#PBS -l nodes=1:ppn=1

module load openmpi/2.1.6
module load scale/6.3.b12

cd $PBS_O_WORKDIR

scalerte -m %%%input_string%%%
grep -a "final result" %%%input_string%%%.out > %%%input_string%%%_done.dat"""

    def move_to_head_node_and_submit_scale_jobs(self, jobs_to_run):
        current_node = socket.gethostname()
        current_directory = os.getcwd()
        print(current_directory)
        ### Make list of all jobs to submit
        for file in jobs_to_run:
            if file.endswith('.sh') == False:
                file = file + ".sh"
            os.system('ssh -tt necluster.ne.utk.edu "cd ' + current_directory + ' && qsub ' + file + '"')

    def run_mt_clutch_job(self, betas, number_of_cases, file_flag="mt_tsunami_", tsunami_template_file = "tsunami_template_default.inp"):
        ### Building template file for this case
        self.build_template_file_for_tsunami(betas,
                                             template_filename=self.template_file)

        ### building inputs, running, combining sdf into a dictionary
        sdf_dict = self.build_run_and_combine_mt_clutch_runs(template_file_string=tsunami_template_file,
                                                             file_flag=file_flag,
                                                             number_to_make=number_of_cases)
        return sdf_dict

    def build_run_and_combine_mt_clutch_runs(self, template_file_string,
                                             file_flag,
                                             number_to_make):

        ### Cleaning up previous _done files
        for file in os.listdir():
            if "_done" in file:
                os.remove(file)
                print("Deleted:", file)

        ### Building Random Number Seed dictionary
        print("Building", number_to_make, "random number seeds")
        target_string = "%%%random_number%%%"
        replacement_dict = collections.OrderedDict()
        replacement_dict[target_string] = []
        for _ in range(number_to_make):
            replacement_dict[target_string].append(str((_) * 10000000000 + 1))
        dictionary_of_replacements = replacement_dict
        self.input_files = []
        for _ in range(number_to_make):
            template_file = open(template_file_string, 'r')
            output_file_string = str(file_flag) + str(_) + ".inp"
            output_file = open(str(output_file_string), 'w')

            for line in template_file:
                for val in dictionary_of_replacements:
                    if val in line:
                        line = line.replace(val, dictionary_of_replacements[val][_])
                output_file.write(line)
            template_file.close()
            output_file.close()
            self.input_files.append(output_file_string)
            print("Built file:", output_file_string)

        ### Building shell script files
        for input_str in self.input_files:
            shell_file_str = input_str + ".sh"
            shell_script = open(shell_file_str, 'w')
            write_string = self.singlethreaded_scale_script.replace("%%%input_string%%%", input_str)
            shell_script.write(write_string)
            shell_script.close()

        ### Submitting scale jobs
        self.move_to_head_node_and_submit_scale_jobs(self.input_files)

        ### Waiting on jobs to be complete
        waiting = True
        while waiting:
            done_file_count = 0
            for file in os.listdir():
                if "_done" in file:
                    done_file_count += 1
            if done_file_count != number_to_make:
                print("Jobs completed:", done_file_count, " waiting on: ", number_to_make - done_file_count)
                print("Waiting 15 seconds")
                time.sleep(15)
            if done_file_count == number_to_make:
                waiting = False
                print("All jobs completed, continuing!")

        print("Combining", len(self.input_files), "sdfs into one dictionary")
        ### Combining sdf for completed files
        combined_sensitivity_dict = self.combine_multiple_sdf_dicts_into_one(self.input_files)

        return combined_sensitivity_dict

    def get_scale_generation_count(self, scale_output):
        ### Adding .out to filename if needed
        if scale_output.endswith('.out') == False:
            temp_str = scale_output.split('.')
            scale_output = temp_str[0] + '.out'

        ### Opening output
        in_data = False
        with open(scale_output) as scale_output_file:
            for line in scale_output_file:
                if "     generation   k-effective" in line:
                    in_data = True
                    continue

                if line.strip() == "":
                    in_data = False
                if in_data:
                    line_split = line.split()
                    if len(line_split) != 8:
                        continue
                    generation = line_split[0]
                    elapsed_time = line_split[7]
        return int(generation)

    def combine_multiple_sdf_dicts_into_one(self, list_of_files):

        list_of_sdf_files = []
        for item in list_of_files:
            if item.endswith('.sdf') == False:
                item_split = item.split(".")
                list_of_sdf_files.append(item_split[0] + ".sdf")
            else:
                list_of_sdf_files.append(item)

        neutrons_per_generation = self.neutrons_per_generation
        skip_generations = self.skip_generations

        ### This part of the function creates a list of all of the sdf files and a list of all
        sdf_dicts = []
        generation_list = []
        for file in list_of_sdf_files:
            print("Parsing sdf file:", file)
            sdf_dicts.append(self.sfh.parse_sdf_file_into_dict(file))
            generation_list.append(self.get_scale_generation_count(file) - skip_generations)
        summed_generation_list = sum(generation_list)

        ### Combines the sdf dicts into a single dict using weighted average and weighted
        ### standard deviation
        combined_dict = collections.OrderedDict()
        for location in sdf_dicts[0]:
            combined_dict[location] = collections.OrderedDict()
            for isotope in sdf_dicts[0][location]:
                combined_dict[location][isotope] = collections.OrderedDict()
                combined_dict[location][isotope]['sensitivity'] = 0.0
                combined_dict[location][isotope]['uncertainty'] = 0.0

                for sdf_count, sdf_dict in enumerate(sdf_dicts):
                    combined_dict[location][isotope]['sensitivity'] += float(
                        sdf_dict[location][isotope]['sensitivity']) * int(
                        generation_list[sdf_count]) / summed_generation_list
                    combined_dict[location][isotope]['uncertainty'] += float(
                        sdf_dict[location][isotope]['uncertainty']) * int(
                        generation_list[sdf_count]) / summed_generation_list

        ### Turns the dict of floats into a dict of strings
        for location in combined_dict:
            for isotope in combined_dict[location]:
                for val in combined_dict[location][isotope]:
                    combined_dict[location][isotope][val] = str(combined_dict[location][isotope][val])

        return combined_dict

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
            material_list.append(self.sfh.combine_material_dicts(material_1, material_2, beta))

        material_string_list = []
        for count, material in enumerate(material_list):
            material_string_list.append(
                self.sfh.build_scale_material_string(material, count + material_count_offset, temperature))

        ### Making list of keys
        flag_list = []
        for x in range(len(material_string_list)):
            flag_list.append(flag.replace(flag_replacement_string, str(x)))

        material_dict = self.sfh.make_data_dict(flag_list, material_string_list)

        self.sfh.create_scale_input_given_target_dict(template_file_string, file_name_flag, material_dict)

    def build_template_file_for_tsunami(self, betas, template_filename,
                                        new_template_filename="tsunami_template_default"):

        material_betas = betas
        materials = self.list_of_material_dictionaries
        self.build_scale_input_from_beta(
            material_betas=material_betas,
            material_1=materials[0],
            material_2=materials[1],
            flag="%material_replace%",
            flag_replacement_string='replace',
            template_file_string=template_filename,
            file_name_flag=new_template_filename)






