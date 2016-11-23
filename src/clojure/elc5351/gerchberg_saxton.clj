(ns elc5351.gerchberg_saxton
  (:use clojure.core
        [elc5351 utils fft image])
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.operators :as op]
            [clojure.core.matrix.linear    :refer [norm]])
  (:import [mikera.matrixx Matrix Matrixx]
           [mikera.vectorz AVector Vector Vectorz]
           [org.jtransforms.fft DoubleFFT_1D DoubleFFT_2D]))

(set! *warn-on-reflection* true)

(m/set-current-implementation :vectorz)

(def ^Matrix source (imread_gray "resources/clinton_gray.png"))
(def ^Matrix target (imread_gray "resources/trump_gray.png"))
;; (def ^Matrix source (let [[rows cols] (m/shape target)
;;                           m ^Matrix (m/mul (Matrix/createRandom rows cols) 255.0)]
;;                       m))

;; (def ^Matrix source2 (let [fft_mat ^Matrix (fft2 target)]
;;                        (complexToReal fft_mat :type :mag)))

(defn ^Matrix gerchberg-saxton [^Matrix source ^Matrix target & {:keys [iterations display]
                                                                 :or {iterations 50}}]
  (let [[rows cols] (m/shape source)
        frame ^Matrix (let [m ^Matrix (Matrix/create rows (bit-shift-left cols 1))]
                        (m/assign! (m/submatrix m [[0 rows] [0    cols]]) (scale_matrix source 255.0))
                        (m/assign! (m/submatrix m [[0 rows] [cols cols]]) (scale_matrix target 255.0))
                        m)
        win (when display
              (imsave frame "output/original.png")
              (imshow frame :title "Source and target"))]
   (loop [i    ^Long   (long 0)
          _X_  ^Matrix (fft2 source :complex false)
          _x_  ^Matrix (Matrix/create rows (bit-shift-left cols 1))
          phase_X nil
          phase_x nil]
     (if (< i iterations)
       (let [X_I   ^Matrix (complexToReal _X_ :type :mag)
             X_P   ^Matrix (complexToReal _X_ :type :phase)
             new_X ^Matrix (composeComplex target X_P :frame _x_)
             x     ^Matrix (ifft2 new_X :complex true :copy false)
             x_P   ^Matrix (complexToReal x :type :phase)
             new_x ^Matrix (composeComplex source x_P :frame _X_)]
         (print ". ") (flush)
         (when win
           (m/assign! (m/submatrix frame [[0 rows] [0    cols]]) (scale_matrix (complexToReal x :type :mag)
                                                                               255.0))
           (m/assign! (m/submatrix frame [[0 rows] [cols cols]]) (scale_matrix X_I 255.0))
           (imupdate win frame :title (format "Iteration %d" i))
           (imsave frame (format "output/iteration_%03d.png" i)))
         (recur (unchecked-inc i) (fft2 new_x :complex true :copy false) _x_ X_P x_P))
       [phase_x phase_X]))))

(gerchberg-saxton source target :iterations 100 :display true)
