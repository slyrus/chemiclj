
# chemiclj

chemiclj is a chemistry library for clojure

# A functional approach to cheminformatics

Clojure encourages a functional programming style in which objects are
immutable. chemiclj is a library for representing various chemical
entities -- molecules, atoms, bons, elements, etc... in a functional,
immutable style. To get a feel for the difference between this
approach and the traditional object-oriented approach to a chemistry
library, consider the task of adding a bond to a molecule. One might
expect that a method, say add-bond, would take a molecule and two
atoms as parameters, and modify the molecule such that it now contains
the new bond. In chemiclj, add-bond returns a completely new molecule,
but one which contains the new bond, and all of the previous
information of the molecule. This approach is used throughout chemiclj
and may take a bit of getting used to.

To see an example, let's use the chemiclj API to build a simple molecule, say, methane:

    (let [c1 (make-atom :c "C1")
          h1 (make-atom :h "H1")
          h2 (make-atom :h "H2")
          h3 (make-atom :h "H3")
          h4 (make-atom :h "H4")]
      (def methane  (add-bond
                     (add-bond
                      (add-bond
                       (add-bond
                        (add-bond (make-molecule #{c1 h1 h2 h3 c2 h4 h5 h6}) c1 h1)
                        c1 h2)
                       c1 h3)
                      c1 h4))))

This may seem tedious at first, but will see how to alleviate this
tedium shortly.

# Dependencies

NOTE: make sure that the dependencies are installed be lein deps, or by
putting the appropriate symlinks in a directory named checkouts.

 shortuct - a graph library for clojure

currently, chemiclj requires (in addition to clojure and clojure-contrib):

* shortcut -- there is no official release (not yet anyway) of a
shortcut jar on clojars, or elsewhere, one MUST INSTALL SHORTCUT BY
HAND AND PLACE A SYMLINK TO IT IN THE checkouts DIRECTORY. This
requires leiningen 1.2 or later.

* fnparse3 -- Joshua Choi's fnparse3 is required for SMILES
parsing. However, the version that exists on github, at the time of
writing this, does not have a project.clj file. I have a forked
repository that can be found at:

`git://github.com/slyrus/fnparse.git`
or
`http://github.com/slyrus/fnparse.git`

note that one needs the develop branch from this repo, so, do:

    git clone git://github.com/slyrus/fnparse.git
    cd fnparse
    git checkout develop

and then make sure that you put a link the fnparse directory in the
chemiclj/checkouts directory.

# Getting Started

Here's a quick example of using chemiclj to read in a compound from a
SMILES string and compute the mass and exact-mass of the molecule:

    user> (require '[chemiclj.core :as chem] '[chemiclj.smiles :as smiles])
    nil
    user> (def vanillin (smiles/read-smiles-string "O=Cc1ccc(O)c(OC)c1"))
    #'user/vanillin
    user> (mass vanillin)
    ; Evaluation aborted.
    user> (chem/mass vanillin)
    152.14732
    user> (chem/exact-mass vanillin)
    152.04736
    user> 

# License

chemiclj is released under a BSD-style license

Cyrus Harmon
ch-github at bobobeach.com
Fri Aug  6 14:21:07 2010

