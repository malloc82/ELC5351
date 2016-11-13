(ns elc5351.fft
  (:use clojure.core
        elc5351.utils)
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.linear :refer [norm]])
  (:import [mikera.matrixx Matrix]
           [org.jtransforms.fft DoubleFFT_2D]))

(set! *warn-on-reflection* true)

(m/set-current-implementation :vectorz)

(defn fft2 [^Matrix m & {:keys [complex copy]
                         :or {copy true}}]
  (let [[rows cols]         (m/shape m)
        complex_mat ^Matrix (if complex (if copy (.clone m) m) (realToComplex m rows cols))]
    (.complexForward (DoubleFFT_2D. rows (if complex (bit-shift-right cols 1) cols))
                     (.asDoubleArray complex_mat))
    complex_mat))


(defn ifft2 [^Matrix m & {:keys [copy scaling]
                          :or {copy true
                               scaling true}}]
  (let [[rows cols] (m/shape m)
        mat_ifft (if copy (.clone m) m)]
    (.complexInverse ^DoubleFFT_2D (DoubleFFT_2D. rows (bit-shift-right cols 1))
                     (.asDoubleArray mat_ifft) true)
    mat_ifft))
