# phenotype.hpoa



It is for creating phenotype.hpoa kg graph using phenotype.hpoa file (http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa)

#Step1 preprocessingt the data
python3 Phenotype_hpoa_file_preprocessing.py -data phenotype.hpoa -output phenotype.hpoa


#Step2 create graph
python3 Phenotype_hpoa_GraphCreation.py -data phenotype.hpoa_12.csv \
            -schema Disease_HPO_Phenotype_Schema.csv -output phenotype.hpoa_12
