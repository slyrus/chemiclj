
# chemiclj

chemiclj is a chemistry library for clojure

# A functional approach to cheminformatics

Clojure encourages a functional programming style in which objects are
immutable. chemiclj is a library for representing various chemical
entities -- molecules, atoms, bons, elements, etc... in a functional
style using immutable objects. To get a feel for the difference
between this approach and the traditional object-oriented approach to
a chemistry library, consider the task of adding a bond to a
molecule. One might expect that a method, say `add-bond`, would take a
molecule and two atoms as parameters, and modify the molecule such
that it now contains the new bond. In chemiclj, add-bond returns a
completely new molecule, but one which contains the new bond, and all
of the previous information of the molecule. This approach is used
throughout chemiclj and may take a bit of getting used to.

To see an example, let's use the chemiclj API to build a simple
molecule, say, methane:

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

Rather than building molecules by hand like this, one usually starts
with a concise representation of a molecule, such as a SMILES
string. For instance, the amino acid valine is represented by the
SMILES string `CC(C)C(C(=O)O)N`. To read a smiles string, use the
function `chemiclj.smiles/read-smiles-string`, as in the following
example:

    user> (require ['chemiclj.smiles :as 'smiles])
    nil
    user> (def alanine (smiles/read-smiles-string "CC(C)C(C(=O)O)N"))
    #'user/alanine

and then to write it back out as a SMILES string:

    user> (smiles/write-smiles-string alanine)
    "CC(C)C(N)C(=O)O"

Writing SMILES strings is still under development, but it should be
fully supported in the near future.

# Dependencies

NOTE: make sure that the dependencies are installed be lein deps, or by
putting the appropriate symlinks in a directory named checkouts.

currently, chemiclj requires (in addition to clojure and clojure-contrib):

* shortcut -- there is no official release (not yet anyway) of a
shortcut jar on clojars, or elsewhere, one MUST INSTALL SHORTCUT BY
HAND AND PLACE A SYMLINK TO IT IN THE checkouts DIRECTORY. This
requires leiningen 1.2 or later.

* fnparse3 -- Joshua Choi's fnparse3 is required for SMILES
parsing. Currently, fnparse3 exists in the develop branch:

    git clone git://github.com/joshua-choi/fnparse.git
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

# License

chemiclj is released under a BSD-style license

Cyrus Harmon
ch-github at bobobeach.com
Fri Aug  6 14:21:07 2010

