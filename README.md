# migadmi

**mig**ration and **admi**xtures: method to estimate parameters in a complex admixture graph with multiple and nested admixture events, when the number of source population could be {1,2,3,...}.

<img src="admixture_parameterization.png" width="400">

For admixture schemes of high complexity, migadmi introduces admixture parameters and makes the symbolic inference of covariance matrix between populations. After that, it optimizes parameters with SLSQP (Sequential Least Squares Programming) method.
Delais of the method is described in the paper [Historical routes for diversification of domesticated chickpea inferred from landrace genomics](https://doi.org/10.1101/2021.01.27.428389).

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

The number of source populations for each mixture could be {1,2,3,4...}. 1 - when it is just a new branch. In principle, the number is not limited; however, the problem could be weakly identifiable with the growth of the number of sources for an admixture event.
To cope with the identification problem, we use regularisation for admixture weights: the Dirichlet prior, with concentration parameter `alpha`, `alpha=1` - the absence of regularisation.

There is an option for step-by-step optimisation of admixture events:

```
admixture_steps = [[0, 1], [2, 3], [4, 5]]
```
In this case the exdixture parameters will be optimised in three steps.

## Run migadmi

```
variables, variance_decomposition, weight_sets = migadmi(tree=tree,
              admixtures=admixtures,
              admixture_steps=admixture_steps,
              dist_matrix=d_mx, alpha=1)
```

To demonstrate the test, please run `example.py`

## Ouput data

#### 1. Optimised parameters
Dictionary of parameters and values. There are three types of parameters: branch lengths, admixture weights, and parts of common variance with sources (alpha parameter on the Figure)
Names of these three types start with the letters "t", "w" and "a", respectively. For example:

```
{'t0': 0.0,
 't1': 0.0,
 't2': 0.729128663694431,
 't3': 0.1202425699727737,
 't4': 0.2139704726613705,
 't5': 0.18213312504471313,
 't6': 0.48002745557810955,
 'a7': 0.23128479141352762,
 'w8': 0.9365519173251352,
 'w9': 0.0,
 'w10': 0.06344808267486482}
 ```


#### 2. Decomposition of variance for mixture populations
For each admixture population, **migadmi** estimates the proportion of variance explained by sources and proportion of own variance.

```
[['ETHI_d', 'LEB_d', 0.0665210082552726],
 ['ETHI_d', 'TUR_d', 0],
 ['ETHI_d', 'IND_d', 8.95945904316835e-5],
 ['ETHI_d', 'ETHI_d', 0.933389397154296]]  # own variance
 ```

#### 3. Names of weight parameters for admixtures
Names of weight parameters, corresponding to the admixture events.

```
[[w8, w9, w10], [w13, w14], [], [w19, w20], [w23, w24], [w27, w28]]
```


## Getting Started

Clone this directory to your computer

## Requirements

To run **migadmi** methods, you need Python 3.4 or later. A list of required Python packages that the migadmi depends on, are in `requirements.txt`.


## Authors

Anna Igolkina developed the **migadmi** package, [e-mail](mailto:igolkinaanna11@gmail.com).

## Citation
This method is described in the paper [Historical routes for diversification of domesticated chickpea inferred from landrace genomics](https://doi.org/10.1101/2021.01.27.428389).


## License information
The **migadmi** package is open-sourced software licensed under the [MIT license](https://opensource.org/licenses/MIT).
