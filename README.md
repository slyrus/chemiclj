
chemiclj is a chemistry library for clojure

NOTE: make sure that the dependencies are installed be lein deps, or by
putting the appropriate symlinks in a directory named checkouts.

Dependencies

 shortuct - a graph library for clojure

currently, chemiclj requires (in addition to clojure and clojure-contrib):

 * shortcut -- there is no official release (not yet anyway) of a
shortcut jar on clojars, or elsewhere, one MUST INSTALL SHORTCUT BY
HAND IN PLACE A SYMLINK TO IT IN THE checkouts DIRECTORY. This
requires leiningen 1.2 or later.

 * fnparse3 -- Joshua Choi's fnparse3 is required for SMILES
 parsing. However, the version that exists on github, at the time of
 writing this, does not have a project.clj file. I have a forked
 repository that can be found at:

   git@github.com:slyrus/fnparse.git

note that one needs the develop branch from this repo, so, do:

    git clone git://github.com/slyrus/fnparse.git
    cd fnparse
    git checkout develop

and then make sure that you put a link the fnparse directory in the
chemiclj/checkouts directory.

Cyrus Harmon
ch-github at bobobeach.com
Fri Aug  6 14:21:07 2010

