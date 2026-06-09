d <- 'dag {

n3 [label="?", fillcolor="#EA3E3E", fontcolor="black"];
n6 [label="Tutor System: Access", fillcolor="#7FD4FF", fontcolor="black"];
n8 [label="Teacher behaviour", fillcolor="#7FD4FF", fontcolor="black"];
n9 [label="Math Knowledge", fillcolor="#A97FFF", fontcolor="black"];
n10 [label="Math Grade", fillcolor="#7FD4FF", fontcolor="black"];
n11 [label="Tutor system: Usage", fillcolor="#7FD4FF", fontcolor="black"];
n12 [label="Curriculum", fillcolor="#7FD4FF", fontcolor="black"];
n13 [label="Student engagement", fillcolor="#A97FFF", fontcolor="black"];
n14 [label="Other maths work", fillcolor="#7FD4FF", fontcolor="black"];
n15 [label="SRL Intervention", fillcolor="#EA9D51", fontcolor="black"];
n16 [label="Other tuition", fillcolor="#7FD4FF", fontcolor="black"];

n9 -> n10 
n6 -> n11 
n11 -> n9 
n6 -> n8 
n8 -> n9 
n8 -> n11 
n13 -> n11 
n13 -> n14 
n14 -> n9 
n15 -> n8 
n15 -> n11 
n11 -> n14 
n11 -> n16 
n16 -> n14 
}
'
# Manually fixing the dag string
d.fixed <- 'dag {

TS.A [label="Tutor System: Access", fillcolor="#7FD4FF", fontcolor="black"];
T.B [label="Teacher behaviour", fillcolor="#7FD4FF", fontcolor="black"];
M.K [label="Math Knowledge", fillcolor="#A97FFF", fontcolor="black"];
M.G [label="Math Grade", fillcolor="#7FD4FF", fontcolor="black"];
TS.U [label="Tutor system: Usage", fillcolor="#7FD4FF", fontcolor="black"];
C [label="Curriculum", fillcolor="#7FD4FF", fontcolor="black"];
S.E [label="Student engagement", fillcolor="#A97FFF", fontcolor="black"];
OMW [label="Other maths work", fillcolor="#7FD4FF", fontcolor="black"];
SRL.I [label="SRL Intervention", fillcolor="#EA9D51", fontcolor="black"];
OT [label="Other tuition", fillcolor="#7FD4FF", fontcolor="black"];

M.K -> M.G 
TS.A -> TS.U 
TS.U -> M.K 
TS.A -> T.B 
T.B -> M.K 
T.B -> TS.U 
S.E -> TS.U 
S.E -> OMW 
OMW -> M.K 
SRL.I -> T.B 
SRL.I -> TS.U 
TS.U -> OMW 
TS.U -> OT 
OT -> OMW 
}
'

# Loading the dag and having a look

library(dagitty)

d <- dagitty(d.fixed)
plot(d)

# Getting the adjustment sets for the total effect of SRL.I on M.G - should be none needed
# This INCLUDES the effect through OT however, which may not be desirable! 
adjustmentSets(d, exposure = "SRL.I", outcome = "M.G")

