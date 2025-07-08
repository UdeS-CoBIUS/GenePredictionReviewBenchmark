# An Overview of Computational Methods for Gene Prediction in Eukaryotes: Strengths, Limitations, and Future Directions. 


<!-- OVERVIEW -->
<h3 id="overview"> Overview </h3>

[To complete. Can add an image if necessary]

> :busts_in_silhouette: __Authors__* `Abigail Djossou and al.`, CoBIUS LAB, Department of Computer Science, Faculty of Science, Université de Sherbrooke, Sherbrooke, Canada*

> :bulb: How to cite us: `upcoming`

> :e-mail: `Contact: abigail.djossou@usherbrooke.ca`

1. [➤ Overview](#overview)
2. [➤ Operating System](#os)
3. [➤ Requirements](#requirements)
4. [➤ Getting Started](#getting-started)
    1. [➤ Example for testing your environment](#main)
    2. [➤ Benchmark](#main)

![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/rainbow.png)

<!-- Operating System -->
<h3 name="os">Operating System</h3>
The benchmarck program was both developed and tested on a system operating Ubuntu version 24.04.6 LTS.

![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/rainbow.png)

<!-- Requirements -->
<h3 id="requirements"> Requirements</h3>

*   __`python3 (at leat python 3.6)`__
*   __`Biopython`__
*   __`bcbio-gff`__

<!-- Getting started -->
<h2 id="getting-started"> :rocket: Getting Started with the benchmark</h2>

* Step 1: Clone the repository with GitHub CLI. or using the command 
<pre><code> git clone https://github.com/UdeS-CoBIUS/GenePredictionReviewBenchmark.git</code></pre>

* Step 2: Install the prerequisites tools with the command <pre><code> ... </code></pre>

* Step 3: Retrieve data from and source codes from repository. You don't have to follow the instructions in this step 3 if you have already cloned the repository in the step 1. 

    * Access and Download the gitlab project of Nicolas Scalzitti and al <a href="https://forge.icube.unistra.fr/n.scalzitti/g3po"> here </a> or clone the repository with the command 

    <pre><code> git clone https://icube-forge.unistra.fr/n.scalzitti/g3po.git </code></pre>

    * Extract the repository `g3po-main.tar.gz`and replace some files that have been updated by using the following commands:

    <pre><code> 
    #description
    mv utils/launch_prediction2.py g3po-main/src/. </code></pre>



<!-- Test -->
<h3 id="test"> :computer: Example for testing your environment</h3>

<!-- Main benchmark -->
<h3 id="main"> :computer: Benchmark</h3>