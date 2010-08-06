
# chemiclj

chemiclj is a chemistry library for clojure

NOTE: make sure that the dependencies are installed be lein deps, or by
putting the appropriate symlinks in a directory named checkouts.

# Dependencies

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

