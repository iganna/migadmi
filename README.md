# migadmi

**mig**ration and **admi**xtures: method to estimate parameters in a complex admixture graph with multiple and nested admixture events, when the number of source population could be {1,2,3,...}.

<img src="admixture_parameterization.png" width="400">

For admixture schemes of high complexity, migadmi introduces admixture parameters and makes the symbolic inference of covariance matrix between populations. After that, it optimizes parameters with SLSQP (Sequential Least Squares Programming) method.

## Input data

#### 1. A binary tree of relationships between base source populations
For some of the population, the relationship should be provided. Then these populations could be used as sources, and after an admixture event, each new mixture population also could play the role of a source. 
We suggest providing the tree in with .nwk format:

```
((LEB_d,TUR_d),(IND_d,UZB_d));
```

To read the tree, use the following code:
```
from ete3 import Tree
file_tree = 'data/tree_init.nwk'
tree = Tree(file_tree)
popset_init, variables_init = get_populations(tree, pop_names)
```



#### 2. Allele frequencies in populations or the distance matrix
Table with populations in columns, SNPs in rows. Allele frequency values are numbers in [0,1]. These values are used to calculate the distance matrix between populations using the compositional data analysis (CoDA).
Of course, if one wants to avoid CoDA, they can provide the pre-calculated distance matrix between populations.

#### 3. Admixture scheme
We suggest providing the scheme in a file and reading it with `read_admixture` function.
The admixture scheme should be organized as follows. Each row is an admixture event. The first word in a row - the name of the mixture population, then sources:

```
ETHI_d: LEB_d, TUR_d, IND_d
MOR_d: LEB_d, TUR_d
TUR_k: TUR_d
UZB_k: UZB_d, TUR_k
MOR_k: MOR_d, TUR_k
LEB_k: LEB_d, TUR_k
```

The number of source populations for each mixture could be {1,2,3,4...}. 1 - when it is just a new branch. In principle, the number is not limited; however, the problem could be weakly identifiable with the growth of the number of sources for an admixture event. To cope with the identification problem, we use regularisation for admixture weights (Dirichlet prior, all concentration parameters = 0.9).

There is an option for step-by-step optimisation of admixture events:

```
admixture_steps = [[0, 1], [2, 3], [4, 5]]
```
In this case the admixture parameters will be optimised in three steps.


## Getting Started

Clone this directory to your computer

## Pipeline (running the test)

To demonstrate the test, please run `pipeline_wnd.py`

## Requirements

To run **migadmi** methods, you need Python 3.4 or later. A list of required Python packages that the migadmi depends on, are in `requirements.txt`.  


## Authors

Anna Igolkina developed the **migadmi** package, [e-mail](mailto:igolkinaanna11@gmail.com).    


## License information
The **migadmi** package is open-sourced software licensed under the [MIT license](https://opensource.org/licenses/MIT).
