# Loading GCMs

library(tidyverse)
library(dagitty)

### Helper functions -----------

to_duration <- function(x){
    x |> 
        str_remove(":00$") |> 
        ms()
}


### Session data ---------------

session_data_raw <- tibble::tribble(
    ~`Session`, ~Stat.familiarity, ~GCM.familiarity, ~GCM.validity, ~Tutorial.duration, ~Modelling.duration,
    "S1 - R1",                 1,               1L,             3,          "4:39:00",          "50:10:00",
    "S1 - R2",                 3,               2L,           3.5,                 NA,                  NA,
    "S1 - R3",                 1,               1L,             3,                 NA,                  NA,
    "S2",                 2,               2L,           3.5,          "4:20:00",          "46:06:00",
    "S3",               2.2,               1L,             3,          "5:53:00",          "36:37:00",
    "S4",                 2,               3L,           1.5,          "4:34:00",          "34:43:00",
    "S5",                 2,               1L,             4,          "2:45:00",          "50:23:00",
    "S6",                 2,               2L,           2.5,          "4:28:00",          "41:13:00",
    "S7",              2.75,               2L,           2.5,          "5:00:00",          "38:55:00",
    "S8",                 4,               2L,             3,          "4:53:00",          "46:04:00"
)

session_data <- tibble::tribble(
    ~`Model`, ~Stat.familiarity, ~GCM.familiarity, ~GCM.validity, ~Tutorial.duration, ~Modelling.duration,
    "S1",                 1.67,               1.33,             3.167,          "4:39:00",          "50:10:00",
    "S2",                 2,               2,           3.5,          "4:20:00",          "46:06:00",
    "S3",               2.2,               1,             3,          "5:53:00",          "36:37:00",
    "S4",                 2,               3,           1.5,          "4:34:00",          "34:43:00",
    "S5",                 2,               1,             4,          "2:45:00",          "50:23:00",
    "S6",                 2,               2,           2.5,          "4:28:00",          "41:13:00",
    "S7",              2.75,               2,           2.5,          "5:00:00",          "38:55:00",
    "S8",                 4,               2,             3,          "4:53:00",          "46:04:00"
) |> 
    mutate(across(ends_with("duration"), to_duration)) 


### Graphical Causal Models --------------

# CD: Course Data (both LMS and Course structure)
# CD.Str Course structure (usually the course data on graphs)
# CD.LMS: LMS activity data
# S. Student attribute. Subsets S.Kn knowledge, S.Cp capabilities, S.Att attribute that is generally unchangible (i.e. gender)
# T. Teacher attribute
# LD. Learning design aspect
# H. Help of some kind. (i.e. assistance)
# Use S.Kn for student knowledge, but S.Kn.Att for a sub component of knowledge
# .Act means "action" or "activity", where people do things or make decisions.

# S.Att.SocCap - social capital
# T.Aff.Integ - teacher (or student) acting with integrity (academic)
# LD.Task - assessment design quality
# T.Cp.Mark - teacher reliable marking
# T/S.Aff.WrkLd - teacher / student workload
dag_s1 <- dagitty('dag S1 {
CD [pos="-0.228,0.457"]
Grade [outcome,pos="-0.830,0.442"]
H [pos="0.285,-0.845"]
LD.Task [pos="-0.408,-0.241"]
S.Aff.Integ [pos="-0.202,-0.798"]
S.Att.SocCap [pos="-0.627,-1.067"]
S.Kn [exposure,pos="0.046,-0.298"]
S.Aff.WrkLd [pos="0.115,0.191"]
T.Aff.Integ [pos="-0.748,-0.410"]
T.Cp.Mark [pos="-1.079,-0.201"]
T.Aff.WrkLd [pos="-1.386,-0.453"]
H -> S.Kn
LD.Task -> CD
LD.Task -> T.Cp.Mark
S.Aff.Integ -> LD.Task
S.Att.SocCap -> H
S.Att.SocCap -> LD.Task
S.Att.SocCap -> T.Aff.Integ
S.Att.SocCap -> T.Cp.Mark
S.Kn -> LD.Task
S.Kn -> S.Aff.WrkLd
S.Aff.WrkLd -> Grade
T.Aff.Integ -> LD.Task
T.Aff.Integ -> T.Cp.Mark
T.Cp.Mark -> CD
T.Cp.Mark -> Grade
T.Aff.WrkLd -> T.Cp.Mark
}')
# plot(dag_s1)

# S.Act - student activity. Doing work
# T.Cp.Inst - teacher instruction (quality of)
# H.Peer - student peer learning / assistance
# T.Fb - teacher feedback
# S.SRL - student SRL

dag_s2 <- dagitty('dag S2 {
CD.LMS [pos="0.190,-0.154"]
CD.Str [pos="-1.216,-0.755"]
Grade [outcome,pos="-0.480,0.370"]
H.AI.evil [pos="-1.007,0.166"]
H.AI.good [pos="-1.043,-0.363"]
H.Peer [pos="-0.300,-0.561"]
H.Tut.Fb [pos="-1.245,-0.054"]
LD [pos="-0.670,-0.848"]
S.Act [pos="-0.608,-0.348"]
S.Cp.SRL [pos="-0.614,0.097"]
S.Kn [exposure,pos="-0.228,0.018"]
T.Cp.Inst [pos="0.073,-0.694"]
H.AI.evil -> Grade
H.AI.good -> CD.Str
H.AI.good -> S.Cp.SRL
H.AI.good -> S.Kn
H.Peer -> S.Kn
H.Tut.Fb -> H.AI.evil
H.Tut.Fb -> S.Cp.SRL
H.Tut.Fb -> S.Kn
LD -> CD.Str
LD -> H.AI.good
LD -> H.Peer
LD -> S.Act
LD -> T.Cp.Inst
S.Act -> CD.LMS
S.Act -> CD.Str
S.Act -> H.Tut.Fb
S.Act -> S.Cp.SRL
S.Act -> S.Kn
S.Cp.SRL -> Grade
S.Cp.SRL -> H.AI.evil
S.Kn -> Grade
S.Kn -> S.Cp.SRL
T.Cp.Inst -> S.Kn
}
')
# plot(dag_s2)

# T.KnSt - Teacher knowing students
# T.AIlit - Teacher AI literacy
# T.Act.UseAI - Teacher AI use
# T.KnSub - Teacher subject knowledge
# S.Act - Student activity, doing the work. Was "S: Course Materials" on graph
# H.Tut - Help from the tutor
# H.SubC - Help from the subject coordinator
# S.Act.Perf - Student performance on the task
# S.Cp.Aca - student academic skill

dag_s3 <- dagitty('dag S3 {
CD.LMS [pos="-0.156,0.899"]
CD.Str [pos="-1.134,0.874"]
Grade [outcome,pos="-0.670,0.870"]
H.AI [pos="0.472,-0.424"]
H.Peer [pos="-0.160,-0.215"]
H.Tut [pos="-0.003,-0.600"]
LD.Task [pos="-1.245,-0.050"]
S.Act [pos="-0.621,-0.571"]
S.Kn [exposure,pos="-0.706,-0.147"]
S.Cp.Aca [pos="0.305,0.083"]
S.Act.Perf [pos="-0.660,0.370"]
T.Act.UseAI [pos="-0.905,-0.032"]
T.Kn.AI [pos="-1.380,0.414"]
T.Kn.Ass [pos="-1.720,0.137"]
T.Kn.Sub [pos="-1.301,-0.830"]
T.Cp.KnStu [pos="-1.638,-0.366"]
H.AI -> H.Tut
H.AI -> S.Act.Perf
H.Peer -> CD.LMS
H.Peer -> S.Kn
H.Peer -> S.Cp.Aca
H.Tut -> CD.LMS
H.Tut -> S.Kn
LD.Task -> CD.Str
LD.Task -> S.Act.Perf
S.Act -> CD.LMS
S.Act -> S.Kn
S.Kn -> S.Act.Perf
S.Cp.Aca -> S.Act.Perf
S.Act.Perf -> Grade
S.Act.Perf -> H.Tut
T.Act.UseAI -> LD.Task
T.Kn.AI -> LD.Task
T.Kn.Ass -> LD.Task
T.Kn.Sub -> H.Tut
T.Kn.Sub -> LD.Task
T.Kn.Sub -> S.Act
T.Cp.KnStu -> LD.Task
T.Cp.KnStu -> S.Kn
}')
# plot(dag_s3)

# S4

# CD.Att - attendance
# S.Cp.Fb - Student feedback literacy
# S.Cp.EvJdg - Student evaluative judgement
# S.Kn.Prior - Base / prior knowledge
# S.Kn - Knowledge to be measured at end of sub
# S.Aff.MindSet - Student affective state: Effort / willingness / learning mindset
# S.Act.Fb - Student acting / taking on feedback
# S.Act.Frm - Student working on formative tasks
# S.Cp.Fb - Student learns from feedback
# S.Att.Nro - Student neuro-divergent 
# S.Att.Hsk - Student attribute, help seeking

dag_s4 <- dagitty('dag S4 {
CD.Att [pos="-1.729,-0.248"]
CD.LMS [pos="-1.344,0.173"]
Grade [outcome,pos="-0.761,0.651"]
H.AI [pos="0.544,-0.025"]
H.Peer [pos="-0.016,-0.496"]
H.Tut [pos="-0.313,-0.694"]
S.Act.Fb [pos="-1.422,-0.686"]
S.Act.Frm [pos="-0.830,-0.686"]
S.Aff.MindSet [pos="-1.671,-1.290"]
S.Att.HSk [pos="0.292,-0.755"]
S.Att.Nro [pos="0.877,-1.078"]
S.Kn [exposure,pos="-0.889,-0.187"]
S.Cp.EvJdg [pos="-0.598,-1.261"]
S.Cp.Fb [pos="-1.026,-1.204"]
S.Kn.Prior [pos="0.014,-1.337"]
S.Cp.Fb [pos="-0.408,-0.172"]
S.Aff.WrkLd [pos="0.436,-1.200"]
H.AI -> Grade
H.AI -> S.Cp.Fb
H.Peer -> S.Cp.Fb
H.Tut -> S.Kn
S.Act.Fb -> S.Act.Frm
S.Act.Frm -> CD.LMS
S.Act.Frm -> S.Kn
S.Aff.MindSet -> CD.Att
S.Aff.MindSet -> S.Act.Fb
S.Aff.MindSet -> S.Act.Frm
S.Att.HSk -> H.Peer
S.Att.HSk -> H.Tut
S.Att.Nro -> H.AI
S.Att.Nro -> S.Att.HSk
S.Kn -> Grade
S.Kn -> S.Act.Fb
S.Cp.EvJdg -> S.Act.Fb
S.Cp.EvJdg -> S.Act.Frm
S.Cp.EvJdg -> S.Cp.Fb
S.Cp.Fb -> S.Act.Fb
S.Cp.Fb -> S.Act.Frm
S.Kn.Prior -> S.Act.Frm
S.Cp.Fb -> S.Kn
S.Aff.WrkLd -> H.AI
S.Aff.WrkLd -> S.Act.Frm
}')
# plot(dag_s4)

# S.Cp.Lrn - Student Learns, in this case used for "knowledge retention"
# S.Aff.Safe - Student feels safe
# S.Aff.Joy - Student is enjoying learning
# S.Aff.MindSet - Student attitude towards learning - change effort above
dag_s5 <- dagitty('dag S5 {
Grade [outcome,pos="-0.480,1.039"]
Grade.Exam [pos="0.105,0.590"]
Grade.Proc [pos="-0.663,0.489"]
Grade.Quiz [pos="-0.323,0.439"]
Grade.Ref [pos="-0.935,0.618"]
H.AI.Agent [pos="-1.308,0.931"]
H.AI.Fb [pos="-1.442,0.618"]
H.AI.Res [pos="-1.252,0.316"]
H.AI.Study [pos="0.472,0.079"]
H.AI.evil [pos="0.563,1.140"]
H.Ment [pos="-0.729,-0.801"]
H.Tut [pos="-0.186,-0.852"]
S.Aff.Joy [latent,pos="-0.130,-0.579"]
S.Aff.MindSet [latent,pos="-0.882,-0.334"]
S.Aff.Safe [latent,pos="-1.111,-0.453"]
S.Kn [exposure,pos="-0.284,-0.223"]
S.Cp.Aca [pos="-0.542,0.054"]
S.Cp.Lrn [exposure,pos="-0.045,0.054"]
Grade.Exam -> Grade
Grade.Proc -> Grade
Grade.Quiz -> Grade
Grade.Ref -> Grade
H.AI.Agent -> Grade.Ref
H.AI.Fb -> Grade.Ref
H.AI.Res -> Grade.Ref
H.AI.Study -> S.Cp.Lrn
H.AI.evil -> Grade.Proc
H.AI.evil -> Grade.Quiz
H.AI.evil -> S.Cp.Lrn
H.Ment -> S.Aff.Joy
H.Ment -> S.Aff.MindSet
H.Ment -> S.Aff.Safe
H.Ment -> S.Cp.Aca
H.Tut -> S.Aff.Joy
H.Tut -> S.Aff.MindSet
H.Tut -> S.Kn
H.Tut -> S.Cp.Aca
S.Aff.Joy -> S.Cp.Lrn
S.Aff.MindSet -> Grade.Ref
S.Aff.MindSet -> S.Kn
S.Aff.MindSet -> S.Cp.Aca
S.Aff.MindSet -> S.Cp.Lrn
S.Aff.Safe -> S.Aff.MindSet
S.Kn -> Grade.Proc
S.Kn -> Grade.Quiz
S.Kn -> S.Cp.Lrn
S.Cp.Aca -> Grade.Exam
S.Cp.Aca -> Grade.Proc
S.Cp.Lrn -> Grade.Exam
}')

# Slightly adjusted to push feedback loops into post framing
# This is used for the effects after the post variable, S.Act.Fb.Post

dag_s6 <- dagitty('dag S6 {
CD.LMS [pos="-0.081,-1.053"]
CD.Str [pos="-0.444,-1.085"]
Grade [outcome,pos="-0.604,0.363"]
H [pos="-0.850,-0.352"]
H.AI.evil [pos="-0.176,-0.305"]
H.AI.good [pos="-0.111,-0.751"]
LD.AIguide [pos="0.223,-0.607"]
LD.Task [pos="-0.817,-0.916"]
S.Act.Fb.Post [pos="-1.366,-0.032"]
S.Act [pos="-0.617,-0.485"]
S.Act.Fb [pos="-0.454,-0.740"]
S.Att.Eq [pos="-1.275,-1.024"]
S.Kn [exposure,pos="-0.853,-0.690"]
S.Cp.Aca [pos="-1.000,-0.039"]
CD.Str -> LD.Task
Grade -> S.Act.Fb.Post
H -> Grade
H -> H.AI.evil
H -> S.Cp.Aca
H.AI.evil -> Grade
H.AI.good -> CD.LMS
H.AI.good <-> S.Act
LD.AIguide -> H.AI.evil
LD.AIguide -> H.AI.good
LD.Task -> S.Kn
LD.Task <-> S.Att.Eq
S.Act.Fb.Post -> S.Kn.Post
S.Act.Fb.Post -> S.Cp.Aca.Post
S.Kn -> S.Kn.Post
S.Cp.Aca -> S.Cp.Aca.Post
S.Act -> Grade
S.Act -> S.Act.Fb
S.Act.Fb -> H.AI.evil
S.Act.Fb -> H.AI.good
S.Act.Fb -> S.Kn
S.Att.Eq -> H
S.Kn -> S.Act
S.Cp.Aca -> Grade
}
')

# Should essentially exclude these two?? They are the post event feedbacks returning to the initial state
# S.Act.Fb.Post -> S.Kn
# S.Act.Fb.Post -> S.Cp.Aca

# S.Kn.Coms - communication skill knowledge / writing / presentation
# S.Kn.Tech - Technical knowledge
# S.Kn.Interp - Interpretation knowledge / skill
# S.Kn.Prior - THis is 'skill' on the original, changed to better match how it is described in other diagrams and this transcript
# S.Act.AIstrat - the student's use choice / strategy for using AI. Higher is a more moral choice. 

dag_s7 <- dagitty('dag S7 {
CD.Att [pos="-0.944,-0.114"]
CD.LMS [pos="-0.933,-0.204"]
CD.LMS.code [pos="-1.142,-0.011"]
CD.LMS.write [pos="-0.934,-0.005"]
Grade [outcome,pos="-1.046,0.161"]
Grade.Coms [pos="-1.045,0.059"]
Grade.Interp [pos="-0.989,0.061"]
Grade.Tech [pos="-1.091,0.057"]
H.AI.Fb [pos="-1.111,-0.190"]
H.AI.Study [pos="-1.124,-0.130"]
H.AI.evil.code [pos="-1.068,-0.011"]
H.AI.evil.write [pos="-1.012,-0.014"]
H.other.evil [pos="-0.969,-0.034"]
S.Act.AIstrat [pos="-1.059,-0.164"]
S.Act.Time [pos="-1.020,-0.165"]
S.Aff.MindSet [pos="-1.038,-0.282"]
S.Kn.Coms [exposure,pos="-1.042,-0.040"]
S.Kn.Interp [exposure,pos="-0.990,-0.068"]
S.Kn.Prior [pos="-0.977,-0.267"]
S.Kn.Tech [exposure,pos="-1.091,-0.066"]
S.Aff.WrkLd [pos="-1.087,-0.276"]
Grade.Coms -> Grade
Grade.Interp -> Grade
Grade.Tech -> Grade
H.AI.Fb -> S.Kn.Coms
H.AI.Fb -> S.Kn.Interp
H.AI.Study -> S.Kn.Tech
H.AI.evil.code -> CD.LMS.code
H.AI.evil.code -> Grade.Tech
H.AI.evil.write -> CD.LMS.write
H.AI.evil.write -> Grade.Coms
H.AI.evil.write -> Grade.Interp
H.other.evil -> CD.LMS.write
H.other.evil -> Grade
S.Act.AIstrat -> H.AI.Fb
S.Act.AIstrat -> H.AI.Study
S.Act.AIstrat -> H.AI.evil.code
S.Act.AIstrat -> H.AI.evil.write
S.Act.Time -> CD.Att
S.Act.Time -> CD.LMS
S.Act.Time -> H.other.evil
S.Act.Time -> S.Act.AIstrat
S.Aff.MindSet -> CD.LMS
S.Aff.MindSet -> H.other.evil
S.Aff.MindSet -> S.Act.AIstrat
S.Aff.MindSet -> S.Act.Time
S.Kn.Coms -> Grade.Coms
S.Kn.Coms <-> S.Kn.Interp
S.Kn.Interp -> Grade.Interp
S.Kn.Prior -> H.other.evil
S.Kn.Prior -> S.Kn.Coms
S.Kn.Prior -> S.Kn.Interp
S.Kn.Prior -> S.Kn.Tech
S.Kn.Tech -> Grade.Tech
S.Kn.Tech -> S.Kn.Interp
S.Aff.WrkLd -> S.Act.AIstrat
S.Aff.WrkLd -> S.Act.Time
}')


# S.SRL - using this, but was referred to as "individual capability" in the original
# S.GrpWrk - students' collaboration ability / capability

dag_s8 <- dagitty('dag S8 {
CD.LMS [pos="-1.617,-1.199"]
Grade [outcome,pos="-0.416,0.274"]
H [pos="0.345,-0.828"]
H.AI.Agent [pos="-0.520,-1.045"]
H.AI.evil [pos="-0.759,-0.884"]
H.Tut [pos="0.213,-1.306"]
LD.Task [pos="-1.134,-0.983"]
LD.Task.Code [pos="-1.867,-0.480"]
LD.Task.Grp [pos="-1.249,-0.224"]
LD.Task.Other [pos="-0.591,-0.662"]
LD.Task.Prac [pos="-1.460,-0.550"]
LD.Task.Write [pos="-1.152,-0.598"]
S.Cp.Grp [exposure,pos="-0.859,-0.142"]
S.Cp.Ind [exposure,pos="-0.405,-0.162"]
S.Kn [exposure,pos="0.166,-0.142"]
CD.LMS -> LD.Task.Code
CD.LMS -> LD.Task.Prac
CD.LMS -> LD.Task.Write
H -> S.Cp.Grp
H -> S.Cp.Ind
H -> S.Kn
H.AI.Agent -> CD.LMS
H.AI.Agent -> H
H.AI.evil -> H
H.AI.evil -> LD.Task
H.Tut -> CD.LMS
H.Tut -> H
LD.Task -> LD.Task.Code
LD.Task -> LD.Task.Prac
LD.Task -> LD.Task.Write
LD.Task.Code -> LD.Task.Grp
LD.Task.Grp -> S.Cp.Grp
LD.Task.Other -> S.Cp.Grp
LD.Task.Other -> S.Cp.Ind
LD.Task.Other -> S.Kn
LD.Task.Prac -> LD.Task.Grp
LD.Task.Write -> LD.Task.Grp
S.Cp.Grp -> Grade
S.Cp.Ind -> Grade
S.Kn -> Grade
}
')

