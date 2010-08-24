
(ns chemiclj.protocol)

(defprotocol PElement
  (isotopes [obj]))

(defprotocol PMass
  (mass [mol])
  (exact-mass [mol]))

(defprotocol PAtomContainer
  (atoms [obj])
  (get-atom [obj id]))

(defprotocol PMolecule
  (add-atom [mol atom])
  (remove-atom [mol atom])
  (bonds [mol] [mol atom])
  (bond? [mol atom1 atom2])
  (add-bond [mol bond] [mol atom1 atom2])
  (remove-bond [mol bond] [mol atom1 atom2])
  (configurations [mol])
  (add-configuration [mol configuration]))


