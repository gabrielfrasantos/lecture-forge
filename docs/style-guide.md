# Style guide

## Teaching pattern
Observation → physical model → derivation → consequence → prediction the viewer can check.

## Figure tags (top right, always visible)
| Tag | Colour | Use |
| --- | --- | --- |
| CALCULATED / SIMULATED | green | Values or motion computed from the stated model; vectors to scale within the figure |
| ANALYTIC | green | Closed-form derivation and its assumptions |
| SCHEMATIC | amber | Illustrative geometry or idealized vectors; never presented as FEA |
| HISTORICAL RECONSTRUCTION | violet | Illustration of an experiment, not the original apparatus |

Field-line arrowheads are static: they show direction, not flow.

## Palette (background #0E1116)
| Meaning | Colour |
| --- | --- |
| Current | #FF9F1C |
| Calculated field | #4CC9F0 |
| Total / resultant field | #BDEFFF (dashed in vector sums) |
| Earth's field | #B388FF |
| Force, torque | #FF4D8D |
| Motion | #7CF29C |
| Moment, space vectors | #FFFFFF (space vectors dashed) |
| Phases a / b / c | #FF6B6B / #FFD166 / #6C8EFF |
| Magnet N / S | #E5484D / #3F7BD9 |

Labels use Inter; equations use LaTeX (Computer Modern). Keep essential content above the bottom 150 px, where players draw captions.

## Voice and audio
- Kokoro-82M, voice `af_heart`, speed 0.98. Write for speech: contractions, connected sentences, equations read as meaning.
- Spoken variables: give single letters letter-phonemes (`[a](/ˈA/)`, `[x](/ˈɛks/)`); say "I times dl", not "I d-l".
- −16 LUFS integrated, true peak below −1 dBTP, no music.
- Captions: at most 2 lines × 42 characters, below 20 characters per second.
