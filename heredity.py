import csv
import itertools
import sys

# python heredity.py data\family0.csv

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])
    print('people:')
    print(people)
    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)

    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                print('have_trait', have_trait)
                print('one_gene', one_gene)
                print('two_genes', two_genes)
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)
                print('- - - - - - - - - -')

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]

def probab(mother_genes, father_genes, child_genes):
    print(f"({mother_genes, father_genes, child_genes})")
    """
    Return probability of number of child genes given
    number of mother and father genes
    """
    # Probabilities that child gets gene from parent
    # parent has 0 genes
    p_from_0_to_1 = 0.01
    # parent has 1 genes
    p_from_1_to_1 = 0.5
    # parent has 2 genes
    p_from_2_to_1 = 0.99

    # Dictionary that map keys (mother_genes, father_genes, child_genes) to probability
    prob = {
        (0,0,0) : 0.99 * 0.99,
        (0,0,1) : 0.01 * 0.99 + 0.99 * 0.01,
        (0,0,2) : 0.01 * 0.01,

        (0,1,0) : 0.99 * 0.5,
        (0,1,1) : 0.01 * 0.5 + 0.99 * 0.5,
        (0,1,2) : 0.01 * 0.5,

        (1,1,0) : 0.5 * 0.5,
        (1,1,1) : 0.5 * 0.5 + 0.5 * 0.5,
        (1,1,2) : 0.5 * 0.5,

        (0,2,0) : 0.99 * 0.01,
        (0,2,1) : 0.99 * 0.99 + 0.01 * 0.01,
        (0,2,2) : 0.01 * 0.99,

        (1,2,0) : 0.5 * 0.01,
        (1,2,1) : 0.5 * 0.01 + 0.5 * 0.99,
        (1,2,2) : 0.5 * 0.99,

        (2,2,0) : 0.01 * 0.01,
        (2,2,1) : 0.99 * 0.01 + 0.01 * 0.99,
        (2,2,2) : 0.99 * 0.99
    }
    key =  (mother_genes, father_genes, child_genes)
    if key not in prob:
        key = (father_genes, mother_genes, child_genes)
    return prob[key]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    print('people dict:')
    print(people)
    names = set(people.keys())
    #print('Names:', names)
    zero_genes = names - one_gene - two_genes
    print('Zero genes', zero_genes)
    print('One gene', one_gene)
    print('Two genes', two_genes)
    n_genes_dict = dict()
    for n in zero_genes:
        n_genes_dict[n] = 0
    for n in one_gene:
        n_genes_dict[n] = 1
    for n in two_genes:
        n_genes_dict[n] = 2
    print('n_genes_dict:')
    print(n_genes_dict)
    no_trait = names - have_trait
    print('Have trait', have_trait)
    print('No trait', no_trait)
    p = 1
    print('Names:')
    i = 0
    for name in list(names):
        i += 1
        print(f"{i}: {name}")
        # How many genes?
        print('Genes:', n_genes_dict[name])
        # Have trait?
        trait = True if name in have_trait else False
        print('Have trait:', trait)
        # Are there parents of name?
        mother = people[name]['mother']
        father = people[name]['father']
        print(f"mother: {mother}, father: {father}")
        #pr_gene = 1
        if not mother:
            print('No parents listed.')
            pr_gene = PROBS["gene"][n_genes_dict[name]]
            print('pr_gene', pr_gene)
        else:
            print('Parents listed.')
            pr_gene = probab(n_genes_dict[mother], n_genes_dict[father], n_genes_dict[name])
            print('pr_gene', pr_gene)
        p *= pr_gene
        #pr_trait = 1
        pr_trait = PROBS['trait'][n_genes_dict[name]][trait]
        print('pr_trait', pr_trait)
        p *= pr_trait

        print('- - - - - -')
    return p

def update(probabilities, one_gene, two_genes, have_trait, p):
    print('update')
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    #print(probabilities)
    names = set(probabilities.keys())
    print('names', names)
    no_genes = names - one_gene - two_genes
    no_trait = names - have_trait
    print('no_genes', no_genes)
    print('one_gene', one_gene)
    print('two_genes', two_genes)
    print('have_trait',  have_trait)
    print('no_trait', no_trait)
    for name in names:
        if name in two_genes:
            probabilities[name]['gene'][2] += p
        elif name in one_gene:
            probabilities[name]['gene'][1] += p
        else:
            probabilities[name]['gene'][0] += p
    for name in names:
        if name in have_trait:
            probabilities[name]['trait'][True] += p
        else:
            probabilities[name]['trait'][False] += p
    print('AFTER update')
    print(probabilities)


def normalize(probabilities):
    print('NORMALIZE')
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for name in probabilities:
        print('name', name)
        gene_sum = 0
        for n_genes in [2,1,0]:
            print('n_genes', n_genes)
            print(probabilities[name]['gene'][n_genes])
            gene_sum += probabilities[name]['gene'][n_genes]
        print('gene_sum', gene_sum)
        for key in probabilities[name]['gene']:
            #print('key', key, 'val', probabilities[name]['gene'][key])
            probabilities[name]['gene'][key] /= gene_sum
        trait_sum = probabilities[name]['trait'][True] + probabilities[name]['trait'][False]
        probabilities[name]['trait'][True] /= trait_sum
        probabilities[name]['trait'][False] /= trait_sum




if __name__ == "__main__":
    main()
