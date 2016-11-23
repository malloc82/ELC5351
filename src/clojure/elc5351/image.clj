(ns elc5351.image
  (:use clojure.core)
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.operators :as op]
            [clojure.java.io :as io])
  (:import[mikera.matrixx.impl  ImmutableMatrix]
          [mikera.matrixx Matrix]
          [mikera.vectorz Vector]
          javax.imageio.ImageIO
          (java.awt.image BufferedImage Raster WritableRaster)
          (javax.swing JFrame JLabel)
          (java.awt Graphics Dimension Color))
  (:gen-class))

(m/set-current-implementation :vectorz)

(defprotocol CA-Image-Protocol
  (imread_gray [s]))

(extend-protocol CA-Image-Protocol
  java.lang.String
  (imread_gray
    ^Matrix
    [path]
    (imread_gray ^BufferedImage (ImageIO/read (io/input-stream path))))
  java.io.File
  (imread_gray
    ^Matrix
    [imfile]
    (imread_gray ^BufferedImage (ImageIO/read (io/input-stream (.getPath imfile)))))
  java.awt.image.BufferedImage
  (imread_gray
    ^Matrix
    [im]
    (let [raster ^Raster        (.getData   im)
          w      ^Integer       (.getWidth  raster)
          h      ^Integer       (.getHeight raster)]
      (case (.getNumBands raster)
        1     (m/reshape (m/matrix :vectorz (.getPixels raster 0 0 w h (double-array (* w h)))) [w h])
        (2 3 4) (let [m ^Matrix (Matrix/create w h)]
                (dotimes [r h]
                  (dotimes [c w]
                    (let [_color ^Color (Color. (.getRGB im c r))]
                      (.set m r c (+ (double (* (.getRed   _color) 0.2989))
                                     (double (* (.getGreen _color) 0.5870))
                                     (double (* (.getBlue  _color) 0.1140)))))))
                m)))))

;; (ns-unmap *ns* 'imshow)
(defmulti imshow (fn [m & _] (class m)))

(defmethod imshow mikera.matrixx.impl.ImmutableMatrix
  [^ImmutableMatrix m & {:keys [scale color title]}]
  (let [_m ^Matrix (m/mutable m)]
    (imshow _m :scale scale :color color :title title)))

(defmethod imshow mikera.matrixx.Matrix
  [^Matrix m & {:keys [scale color title]}]
  (let [[^Integer row ^Integer col] (m/shape m)
        _m    ^Matrix        (m/clone m)
        image ^BufferedImage (if color
                               (BufferedImage. col row BufferedImage/TYPE_INT_ARGB)
                               (BufferedImage. col row BufferedImage/TYPE_BYTE_GRAY))]
    ;; (println "in imshow matrix")
    (when scale
      (let [_max ^Double (.elementMax _m)
            _min ^Double (.elementMin _m)
            _r   ^Double (- _max _min)]
        (if (not= _r 0.0)
         (m/emap! #(Math/round ^Double (double (* (/ (- % _min) _r) scale))) _m))))
    (if color
      (let [len ^Integer  (* row col)
            arr  ^doubles (m/to-double-array _m nil)
            argb ^ints    (int-array len)]
        (dotimes [i len]
          (aset argb i (Color/HSBtoRGB (/ (- 255.0 (aget arr i)) 360.0) 1.0 1.0)))
        (.setRGB image 0 0 col row argb 0 row))
      (.setPixels ^WritableRaster (.getRaster image) 0 0 col row (m/to-double-array _m nil)))
    (imshow image :title title)))

(defmethod imshow java.awt.image.BufferedImage
  [^BufferedImage im & {:keys [title]}]
  ;; (println "in imshow bufferedimage")
  (let [canvas (proxy [JLabel] []
                 (paint [g] (.drawImage ^Graphics g im 0 0 this)))]
    {:jframe (doto (JFrame.)
               (.add canvas)
               (.setSize (Dimension. (.getWidth  im) (+ 22 (.getHeight im))))
               ;; (.setSize (Dimension. (+ (.getWidth  im) 15) (+ (.getHeight im) 30)))
               (.setTitle ^String (if title (str title) "no name"))
               (.show)
               (.toFront)
               (.setVisible true))
     :image im}))

(defn imupdate
  [{^JFrame frame :jframe ^BufferedImage image :image} ^Matrix m & {:keys [scale title]}]
  (let [row ^Integer (.getHeight image)
        col ^Integer (.getWidth  image)
        ;; [^Integer row ^Integer col] (m/shape m)
        _m ^Matrix (m/clone m)
        s  ^Double (.elementMax _m)]
    (when scale (m/emap! #(if (= % 0.0) 0.0 (Math/round ^Double (double (* % (/ scale s))))) _m))
    (if (or (== (.getType image) BufferedImage/TYPE_INT_ARGB)
            (== (.getType image) BufferedImage/TYPE_INT_RGB))
      (let [len ^Integer  (* row col)
            arr  ^doubles (m/to-double-array _m nil)
            argb ^ints    (int-array (* len))]
        (dotimes [i len]
          (aset argb i (Color/HSBtoRGB (/ (- 255.0 (aget arr i)) 360.0) 1.0 1.0)))
        (.setRGB image 0 0 col row argb 0 row))
      (.setPixels ^WritableRaster (.getRaster image) 0 0 col row ^doubles (m/to-double-array _m nil)))
    (.repaint frame)
    (when title (.setTitle frame (str title)))
    {:jframe frame
     :image image}))

(defn imsave [^Matrix im ^String filename & {:keys [scale color ftype]
                                             :or   {color false
                                                    ;; imtype BufferedImage/TYPE_BYTE_GRAY
                                                    ftype  "png"}}]
  (let [[^Integer row ^Integer col] (if (m/matrix? im)
                                      (m/shape im)
                                      (m/shape (im 0)))
        imtype ^Integer        (if color BufferedImage/TYPE_INT_ARGB BufferedImage/TYPE_BYTE_GRAY)
        image  ^BufferedImage  (BufferedImage. col row imtype)
        raster ^WritableRaster (.getRaster image)
        _m     ^Matrix         (m/clone im)]
    (when scale
      (let [_max ^Double (.elementMax _m)
            _min ^Double (.elementMin _m)
            _r   ^Double (- _max _min)]
        (if (not= _r 0.0)
         (m/emap! #(Math/round ^Double (double (* (/ (- % _min) _r) scale))) _m))))
    (if color
      (let [len ^Integer  (* row col)
            arr  ^doubles (m/to-double-array _m nil)
            argb ^ints    (int-array len)]
        (dotimes [i len]
          (aset argb i (Color/HSBtoRGB (/ (- 255.0 (aget arr i)) 360.0) 1.0 1.0)))
        (.setRGB image 0 0 col row argb 0 row))
      (.setPixels ^WritableRaster (.getRaster image) 0 0 col row (m/to-double-array _m nil)))
    (ImageIO/write image ^String ftype (io/file filename))))

(defn im2bw!
  ^Matrix
  [^Matrix m func]
  (let [[^Integer row ^Integer col] (m/shape m)]
    (dotimes [r row]
      (dotimes [c col]
        (.set m r c (func (.get m r c)))))
    m))
