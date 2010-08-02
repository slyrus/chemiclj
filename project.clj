
(defproject chemiclj "0.0.0"
  :description "A computational chemistry library in clojure"
  :dependencies [[org.clojure/clojure "1.2.0-RC1"]
                 [org.clojure/clojure-contrib "1.2.0-RC1"]]
  :dev-dependencies [[swank-clojure "1.3.0-SNAPSHOT"]]
  :jvm-opts ["-Xmx1g -XX:+UseConcMarkSweepGC -XX:+CMSClassUnloadingEnabled -XX:MaxPermSize=256m"])

