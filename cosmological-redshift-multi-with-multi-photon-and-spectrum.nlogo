breed [MWs MW]
breed [galaxies galaxy]
breed [photons photon]

globals [rate-per-tick
  gy-per-tick
  gly-per-pixel
  current-separation
  current-time
  us
  Hubble-plot-data
  should-stop?
  h
  c
  kb
  n-wavelengths
  spectrum
  spectrum0
  wavelengths
  intensities
  wavelength-min
  wavelength-peak
  wavelength-max
  initial-distance
  distance-when-received
]

turtles-own [
  separation
  recession-rate
  my-id
  generation
]

; galaxies-own [
; ]

photons-own [
  made-by-galaxy
  made-at-time
  made-at-distance
  redshift
  wavelength
  wavelength0
  intensity
  CMB?
]

to setup
  clear-all
  random-seed new-seed

  ; constants
  set h 6.626E-34
  set c 2.997E8
  set kb 1.381E-23

  set wavelength-peak 2.8977E6 / transparency-temperature
  set wavelength-min wavelength-peak * 0.5
  set wavelength-max wavelength-peak * 2.57
  set n-wavelengths 2 * max-pycor + 1


  set should-stop? false

  set gy-per-tick 0.1     ; each tick of the clock is 100 million years
  set gly-per-pixel 0.1   ; 10 pixels is 1 billion light-years

  set current-time 0
  set rate-per-tick ln (1 + expansion / 100) * gy-per-tick
  set-default-shape photons "circle"

  ; from galaxy, at time, at distance, initial wavelength, observed wavelength,
  ;   received-time, galaxy-distance at time received
  set Hubble-plot-data (list (list 0 0 0 0 0 0 0))

  set-current-plot "Spectrum from primordial gas"
  set-current-plot-pen "received"
  set-plot-pen-mode 0
  set-current-plot-pen "initial"
  set-plot-pen-mode 0
  set-plot-y-range 0.1 1.05

  reset-ticks

  create-MWS 1
  [                                 ; create one turtle
    set color blue
    set xcor (-0.9) * max-pxcor      ; set the turtle's initial position and heading
    set ycor 0
    set heading 0
    set shape "circle"
    set size 5                      ; easier to see this way
    set us who
    ;pen-down
  ]

  create-galaxies n-galaxies
  [                                 ; create one Milky Way
    set color red
    ;set xcor (-0.87) * max-pxcor + who * 1.87 * max-pxcor / (n-galaxies + 1)   ; set the turtle's initial position and heading
    set xcor (-0.9) * max-pxcor + who / gly-per-pixel   ; set the turtle's initial position and heading
    ;set ycor random-ycor
    set ycor random-ycor * 0.8
    set heading 90
    set shape "circle"
    set size 3                      ; easier to see this way
    set my-id who
    set label who
    set separation (xcor - [xcor] of turtle us) * gly-per-pixel
    set generation 0
  ]

  foreach sort galaxies [x -> ask x
    [
      hatch-photons n-wavelengths [
        set CMB? true
        set label ""
        set heading 270
        set size 1
        set made-by-galaxy my-id
        set made-at-time ticks
        set made-at-distance separation
        set my-id ([who] of self - (n-galaxies + 1) - (made-by-galaxy - 1) * n-wavelengths)
        set ycor (max-pycor - my-id)
        set wavelength wavelength-min + my-id * (wavelength-max - wavelength-min) / (2 * max-pycor + 1)
        set wavelength0 wavelength
        set color wavelength-to-color wavelength
        set intensity planck-function transparency-temperature wavelength
        if (([ycor] of galaxy made-by-galaxy - ycor) > 1) and (([ycor] of galaxy made-by-galaxy - ycor) <= 2) [
          set label made-by-galaxy
          ]
        ]

      set spectrum [intensity] of photons with [made-by-galaxy = [who] of myself]
      ask photons with [made-by-galaxy = [who] of myself] [
        set intensity (intensity / (max spectrum))
        set color scale-color color intensity -0.1 1.5
      ]

      if who = 1 [
        set-current-plot "Spectrum from primordial gas"
        set-current-plot-pen "initial"
        set-plot-pen-mode 0
        set spectrum0 [(list wavelength intensity)] of photons with [made-by-galaxy = [who] of myself]
        set spectrum0 sort-by [ [row1 row2] -> (item 0 row1) < (item 0 row2) ] spectrum0
        ;show spectrum
        set wavelengths map [ row -> item 0 row ] spectrum0
        set intensities map [ row -> item 1 row ] spectrum0
        if output-spectra? [
          show wavelengths
          show intensities
        ]
        (foreach wavelengths intensities
          [ [a b] -> plotxy a b])
      ]
;      hatch-photons 1 [
;        set CMB? false
;        set ycor ycor - 2.5
;        set xcor xcor + 0.5
;        set label my-id
;        set heading 270
;        set size 0
;      ]
    ]  ; end commands for each galaxy doing setup
  ]    ; end loop through all galaxies
  tick
end

to go

  if should-stop? [
    ;show "stop!"
    set-current-plot "Spectrum from primordial gas"
    set-current-plot-pen "received"
    ask photons with [generation = 0 and (xcor <= (-0.9) * max-pxcor - 1) and CMB? and my-id = 0] [
      set initial-distance made-at-distance
      set distance-when-received [separation] of galaxy made-by-galaxy
    ]
    set spectrum [(list wavelength intensity)] of photons with [generation = 0 and (xcor <= (-0.9) * max-pxcor - 1) and CMB?]
    set spectrum sort-by [ [row1 row2] -> (item 0 row1) < (item 0 row2) ] spectrum
    ;show spectrum
    set wavelengths map [ row -> item 0 row ] spectrum
    set intensities map [ row -> item 1 row ] spectrum
    if output-spectra? [
      show wavelengths
      show intensities
    ]
    plot-pen-up
    plotxy item 0 wavelengths item 0 intensities
    plot-pen-down
    (foreach wavelengths intensities
      [ [a b] -> plotxy a b])
    stop
  ]
  update
  ask photons
  [
    if xcor > (-0.9) * max-pxcor [
      set separation (xcor - [xcor] of turtle us) * gly-per-pixel
      set recession-rate rate-per-tick * separation / gly-per-pixel
      fd 1 - recession-rate
      set wavelength (wavelength * (1 + rate-per-tick))
      set color wavelength-to-color wavelength
      set color scale-color color intensity -0.1 1.5
    ]
    if (xcor <= (-0.9) * max-pxcor - 1) and CMB? [
      set should-stop? true
    ]
    if (xcor <= (-0.9) * max-pxcor - 1) and not CMB? [
      die
    ]
    if xcor <= (-0.9) * max-pxcor [
      ; record data to list
      set Hubble-plot-data lput (list made-by-galaxy made-at-time made-at-distance wavelength0 wavelength (ticks * gy-per-tick) ([separation] of turtle made-by-galaxy)) Hubble-plot-data ;
;      plotxy wavelength intensity
      fd 1 - recession-rate
    ]
    if (xcor >= (max-pxcor - 1)) [
      die
    ]
  ]
  ask galaxies
  [
    ;set separation (xcor - [xcor] of turtle us) * gly-per-pixel
    set recession-rate rate-per-tick * separation / gly-per-pixel
    set separation separation + recession-rate * gly-per-pixel
    fd recession-rate
    if xcor >= (max-pxcor - 1) [hide-turtle]
    if ((ticks mod photons-every) = 0) and (xcor < (max-pxcor - 1)) and (not hidden?) and later-photons? [
      set generation generation + 1
      hatch-photons 1 [
        set CMB? false
        set label ""
        set heading 270
        set size 1.5
        set made-by-galaxy my-id
        set made-at-time ticks
        set made-at-distance separation
        set wavelength 400
        set wavelength0 wavelength
        set color wavelength-to-color wavelength
        set intensity 0.7
      ]
      hatch-photons 1 [
        set CMB? false
        set label precision current-time 3
        set heading 270
        set size 0
        set ycor ycor - 2.5
        set made-by-galaxy my-id
        set made-at-time ticks
        set made-at-distance separation
        set wavelength 0
        set color wavelength-to-color wavelength
        set intensity 0
      ]
    ]
  ]
  tick

end

to update
  set rate-per-tick ln (1 + expansion / 100) * gy-per-tick
  set current-time ticks * gy-per-tick
  ;if ticks mod 30 = 0 [reset-my-plots]
end

to-report wavelength-to-color [lambda]
  ; original routine from Andrew Duffy, modified
  let R 0
  if (lambda < 400) [set R 76]
  if ((lambda >= 400) and (lambda <= 500)) [set R floor (160 - 160 * (lambda - 400) / 100)]
  if ((lambda >= 558) and (lambda < 590)) [set R floor (255 - 255 * (590 - lambda) / (32 ^ 2) )]
  if ((lambda >= 590) and (lambda < 650)) [set R 255]
  if ((lambda >= 650) and (lambda <= 700)) [set R floor (255 - 2 * (lambda - 650))]
  if (lambda > 700) [set R 155]
  let G 0
  if (lambda < 460) [set G 49]
  if (lambda >= 460) and (lambda <= 500) [set G floor (255 - 255 * (500 - lambda) ^ 2 / 1600)]
  if (lambda >= 500) and (lambda <= 570) [set G 255]
  if (lambda > 570) and (lambda <= 640) [set G floor (255 - 255 * (lambda - 570) ^ 2 / (70 ^ 2))]
  if (lambda > 640) [set G 0]
  let B 0
  if (lambda < 400) [set B 101]
  if (lambda >= 400) and (lambda < 460) [set B 255]
  if (lambda >= 460) and (lambda < 550) [set B floor (255 - 255 * (lambda - 460) ^ 2 / (90 ^ 2))]
  report approximate-rgb R G B
end



to-report planck-function [temperature lambda]
  ; for wavelength in nm
  let lam lambda * 1E-9
  let a (2 * h * c ^ 2)
  let b lam ^ 5
  let f h * c
  let g lam * kb * temperature
  report ((a / b) / (e ^ (f / g) - 1))

end

to continue
  set should-stop? false
  ask photons [if xcor <= (-0.9) * max-pxcor - 1 [die]]
end

to clean-plot
  set-current-plot "Spectrum from primordial gas"
;  clear-plot
  set-current-plot-pen "received"
  plot-pen-reset
  set-plot-pen-mode 0
;  set wavelengths map [ row -> item 0 row ] spectrum0
;  set intensities map [ row -> item 1 row ] spectrum0
;  (foreach wavelengths intensities
;    [ [a b] -> plotxy a b])
end


; Model developed by Geoff Mathews, as part of the Foothill AstroSims project
@#$#@#$#@
GRAPHICS-WINDOW
0
55
1113
344
-1
-1
5.5
1
10
1
1
1
0
0
0
1
-100
100
-25
25
1
1
1
ticks
30.0

BUTTON
10
10
77
44
Setup
setup
NIL
1
T
OBSERVER
NIL
S
NIL
NIL
1

BUTTON
761
10
867
43
Go
go
T
1
T
OBSERVER
NIL
G
NIL
NIL
0

SLIDER
496
10
707
43
expansion
expansion
0
10
7.0
0.1
1
% per billion years
HORIZONTAL

MONITOR
74
386
344
443
time (billions of years)
current-time
2
1
14

BUTTON
194
702
260
736
Step
go
NIL
1
T
OBSERVER
NIL
O
NIL
NIL
1

TEXTBOX
618
686
1144
736
Foothill AstroSims, https://foothill.edu/astronomy/astrosims.html, CC BY-NC
13
0.0
1

BUTTON
505
703
571
736
print data
print \"from galaxy, at time, at distance, initial wavelength, observed wavelength, received-time, galaxy-distance at time received\"\nprint Hubble-plot-data
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
12
724
186
757
photons-every
photons-every
1
10
5.0
1
1
tick
HORIZONTAL

MONITOR
71
557
343
614
distance when light first received
distance-when-received
2
1
14

SLIDER
271
10
488
43
transparency-temperature
transparency-temperature
1000
10000
3000.0
100
1
K
HORIZONTAL

PLOT
362
349
1114
655
Spectrum from primordial gas
wavelength (nm
relative brightness
400.0
700.0
0.1
1.05
true
true
"set-plot-pen-mode 0\nset-plot-y-range 0 1.05\nset-plot-x-range 400 700" ""
PENS
"initial" 1.0 0 -16777216 true "" ""
"received" 1.0 0 -2674135 true "" ""

SWITCH
12
687
186
720
later-photons?
later-photons?
1
1
-1000

SLIDER
92
10
264
43
n-galaxies
n-galaxies
1
19
18.0
1
1
NIL
HORIZONTAL

BUTTON
871
10
957
43
continue
continue
NIL
1
T
OBSERVER
NIL
C
NIL
NIL
1

BUTTON
1025
10
1114
43
clear plot
clean-plot
NIL
1
T
OBSERVER
NIL
P
NIL
NIL
1

MONITOR
72
490
342
547
Initial distance
initial-distance
2
1
14

SWITCH
344
704
501
737
output-spectra?
output-spectra?
1
1
-1000

@#$#@#$#@
## WHAT IS IT?

This shows all of the code needed to set up a plot which will display two variables as the model runs.  In this case, the two variables are the x-coordinate and the y-coordinate of a turtle.

<!-- 2004 -->
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.4.0
@#$#@#$#@
need-to-manually-make-preview-for-this-model
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
