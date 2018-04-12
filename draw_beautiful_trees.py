#!/usr/bin/env python

from ete3 import Tree,TreeStyle,TextFace

#t = Tree("(((((FelCat:1.0, (PriBen:1.0, PriViv:1.0)Anc2:1.0)Anc3:1.0, (AciJub:1.0, PumCon:1.0)Anc1:1.0)Anc4:1.0, LynPar:1.0)Anc5:1.0, CarCar:1.0)Anc6:1.0, (PanTig:1.0,(PanOnc:1.0, (PanLeo:1.0, PanPar:1.0)Anc7:1.0)Anc8:1.0)Anc9:1.0);", format=1)

#Topolgy from Evolution of Cats, 2007
#t = Tree("(((((FelCat:2.0, (PriBen:1.0, PriViv:1.0)Anc2:1.0), (AciJub:2.0, PumCon:2.0)Anc1:1.0), LynPar:4.0), CarCar:5.0), (PanTig:5.0,(PanOnc:4.0, (PanLeo:3.0, PanPar:3.0)Anc7:1.0)Anc8:1.0)Anc9:1.0);", format=1)

#Topology from Kliver
#t = Tree("(((((FelCat:2.0, (PriBen:1.0, PriViv:1.0)Anc2:1.0), LynPar:3.0), (AciJub:3.0, PumCon:3.0)Anc1:1.0), CarCar:5.0), (PanTig:5.0,(PanOnc:4.0, (PanLeo:3.0, PanPar:3.0)Anc7:1.0)Anc8:1.0)Anc9:1.0);", format=1)

#branch length = Kliver's tree branches lengths*1000
#t = Tree("((((AciJub:2.91010233990525967,PumCon:2.18060580721677255)Anc1:1.05347870516064369,(((PriBen:0.98094949038301073,PriViv:1.00313958496688689)Anc2:1.80011684448748259,FelCat:2.86356101725523731)100:3.3638608143336113,LynPar:3.07709301059199899)96:0.14292306680570912)100:0.59099766777017380,CarCar:3.15189615063413749)100:0.78860645297528301,((PanOnc:1.05156689588516367,(PanLeo:1.24626526455903818,PanPar:0.96868719506379029)Anc7:0.28662050599669259)Anc8:0.38474422414646876,PanTig:1.69403903125740969)Anc9:2.67796389385434753);", format=1)
t = Tree("((CroCro:0.03163103687492756222,((((AciJub:0.00291010233990525967,PumCon:0.00218060580721677255)Anc1:0.00105347870516064369,(((PriBen:0.00098094949038301073,PriViv:0.00100313958496688689)Anc2:0.00180011684448748259,FelCat:0.00286356101725523731)100:0.00033638608143336113,LynPar:0.00307709301059199899)96:0.00014292306680570912)100:0.00059099766777017380,CarCar:0.00315189615063413749)100:0.00078860645297528301,((PanOnc:0.00105156689588516367,(PanLeo:0.00124626526455903818,PanPar:0.00096868719506379029)Anc7:0.00028662050599669259)Anc8:0.00038474422414646876,PanTig:0.00169403903125740969)Anc9:0.00267796389385434753)100:0.01757780609900925356):0.03362414544383966059,CanFam:0.03362414544383966059);", format=1)
t.convert_to_ultrametric()

print t.write(format=1) 


(t & "AciJub").add_features(label="150")
(t & "AciJub").add_face(TextFace((t & "AciJub").label, ftype='Arial'), column=0, position="branch-top")
(t & "PumCon").add_features(label="293")
(t & "PumCon").add_face(TextFace((t & "PumCon").label, ftype='Arial'), column=0, position="branch-top")
(t & "CarCar").add_features(label="185")
(t & "CarCar").add_face(TextFace((t & "CarCar").label, ftype='Arial'), column=0, position="branch-top")
(t & "LynPar").add_features(label="170")
(t & "LynPar").add_face(TextFace((t & "LynPar").label, ftype='Arial'), column=0, position="branch-top")
(t & "PanPar").add_features(label="108")
(t & "PanPar").add_face(TextFace((t & "PanPar").label, ftype='Arial'), column=0, position="branch-top")
(t & "PanLeo").add_features(label="191")
(t & "PanLeo").add_face(TextFace((t & "PanLeo").label, ftype='Arial'), column=0, position="branch-top")
(t & "PanOnc").add_features(label="112")
(t & "PanOnc").add_face(TextFace((t & "PanOnc").label, ftype='Arial'), column=0, position="branch-top")
(t & "PanTig").add_features(label="181")
(t & "PanTig").add_face(TextFace((t & "PanTig").label, ftype='Arial'), column=0, position="branch-top")
(t & "PriViv").add_features(label="137")
(t & "PriViv").add_face(TextFace((t & "PriViv").label, ftype='Arial'), column=0, position="branch-top")
(t & "PriBen").add_features(label="392")
(t & "PriBen").add_face(TextFace((t & "PriBen").label, ftype='Arial'), column=0, position="branch-top")
#PanLeo+PanPar
(t & "Anc7").add_features(label="19 ")
(t & "Anc7").add_face(TextFace((t & "Anc7").label, ftype='Arial'), column=0, position="branch-top")
#PanLeo+PanOnc+PanPar
(t & "Anc8").add_features(label="2 ")
(t & "Anc8").add_face(TextFace((t & "Anc8").label, ftype='Arial'), column=0, position="branch-top")
#PanLeo+PanOnc+PanPar+PanTig
(t & "Anc9").add_features(label="3 ")
(t & "Anc9").add_face(TextFace((t & "Anc9").label, ftype='Arial'), column=0, position="branch-top")
#AciJub+PumCon
(t & "Anc1").add_features(label="5 ")
(t & "Anc1").add_face(TextFace((t & "Anc1").label, ftype='Arial'), column=0, position="branch-top")
#PriBen+PriViv
(t & "Anc2").add_features(label="15 ")
(t & "Anc2").add_face(TextFace((t & "Anc2").label, ftype='Arial'), column=0, position="branch-top")

for node in t.traverse():
    node.img_style['size'] = 0
    if node.is_leaf():
       name_face = TextFace(' '+node.name, ftype='Arial', fstyle='italic', fsize=10)
       node.add_face(name_face, column=0, position="branch-right")

#t.get_ascii(attributes=["name", "label"])
ts = TreeStyle()
ts.scale = 12000 # 12000 pixels per branch length unit
#t.show(tree_style=ts)
ts.show_leaf_name = False
ts.branch_vertical_margin = 10
#ts.show_branch_support = True
#t.get_ascii(attributes=["name", "label"])
t.render('tree.pdf', tree_style=ts)
