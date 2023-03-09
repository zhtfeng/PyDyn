#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 19:52:06 2019

@author: zhtfeng
"""


class Configuration:
    # all configurations were stored in a dict that you can call with self.config[configname]
    # functions are used in main to read config and output config, if no config is updated, a default config will
    # be used
    def __init__(self, config_file, output_file, output_level):

        self.file = config_file
        self.config = { \
 \
            "charge": 0, \
            "multiplicity": 1, \
            "memory": 60, \
            "processors": 10, \
            "chk": " ", \
            "classical": 1, \
            "numimage": 0, \
            "method": "B3LYP/6-31G*", \
            "title": "Untitled", \
            "cannonball": False, \
            "twodirection": False, \
            "diagnostic": False, \
            "timestep": 1., \
            "scaling": 1., \
            "temperature": 273.15, \
            "Initialdis": 1,
            "classicalSpacing": 2}

        self.output = str(
            output_file)  # set up the name of the output file, the file will be created in the print function
        self.output_level = output_level

    def read_config(self):  # read the config file and update corresponding setup into self.config

        open_file = open(self.file)
        config_file_lines = open_file.readlines()
        config_file_lines = [lines for lines in config_file_lines if ":" in lines]
        config_from_read = {}
        for i, each_line in enumerate(config_file_lines):
            
            config_file_lines[i] = each_line.split(':', 1)
            config_file_lines[i][1] = config_file_lines[i][1].title()
            config_file_lines[i][0] = config_file_lines[i][0].strip()

        for each_line in config_file_lines:

            try:

                value = each_line[1]
                config_from_read[str(each_line[0])] = eval(value)

            except SyntaxError:

                config_from_read[str(each_line[0])] = each_line[1].strip()

            except NameError:

                config_from_read[str(each_line[0])] = each_line[1].strip()

        self.config.update(config_from_read)

    def print_config(self):  # call this function for printing out the config info into out.txt

        output_file = open('out.txt', 'w')
        output_file.writelines(['***************** Reading Configuration *****************\n', '\n'])

        for key, value in self.config.items():
            output_file.writelines([str(key), ':        ', str(value), '\n'])

        output_file.writelines(['\n', '***************** End Reading Configuration *****************\n'])
        output_file.close()


test = Configuration("config.txt", 'out.txt', None)
test.read_config()
test.print_config()
