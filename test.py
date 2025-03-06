from heredity import joint_probability, probab
people = {'Lily': {"mother": None, "father": None},
          'James': {"mother": None, "father": None},
          'Harry': {"mother": "Lily", "father": "James"}}
one_gene = {'Harry'}
two_genes = {'James'}
have_trait = {'James'}
p = joint_probability(people, one_gene, two_genes, have_trait)
print(p)
# Result should be = 0.0026643247488


#key = (0,1,1)
#print(probab(*key))
#key = (1,0,1)
#print(probab(*key))