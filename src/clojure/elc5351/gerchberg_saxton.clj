(ns elc5351.gerchberg_saxton
  (:use clojure.core
        elc5351.fft)
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.operators :as op]
            [clojure.core.matrix.linear    :refer [norm]])
  (:import [mikera.matrixx Matrix Matrixx]
           [mikera.vectorz AVector Vector Vectorz]
           [org.jtransforms.fft DoubleFFT_1D DoubleFFT_2D]))

(set! *warn-on-reflection* true)

(m/set-current-implementation :vectorz)
