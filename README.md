# crisprtools

crisprtools is a package to allow for the analysis of CRISPR Screens.
It includes 3 algorithms to analyze CRISPR screens developed by other individuals:

1. drugZ - was developed in the Traver Hart lab, by Gang Wang 
(GitHub Source Code: https://github.com/hart-lab/drugz, Publication: https://www.biorxiv.org/content/early/2017/12/12/232736)

2. BAGEL - was developed in the Traver Hart lab, by Traver Hart 
(GitHub Source Code: https://github.com/hart-lab/bagel, Publication: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1015-8)

3. JACKS - was developed in the Leopold Parts lab, by Fellicity Allen 
(GitHub Source Code: https://github.com/felicityallen/JACKS/tree/master/2018_paper_materials, Publication: https://www.biorxiv.org/content/early/2018/04/12/285114)

It also includes an algorithm to filter out sgRNA's based on their consistency between replicates and conditions. It uses the Bhattacharyya Coefficient to calculate the overlap between two distributions of replicate samples. The bhatt.coeff.R script was written by Thomas Guillerme (GitHub: https://github.com/TGuillerme/Total_Evidence_Method-Missing_data/blob/master/Functions/TreeComparisons/bhatt.coef.R)

Will add more useful information soon...

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

install.packages("devtools")
library(devtools)
install_github("singjc/crisprtools")

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

