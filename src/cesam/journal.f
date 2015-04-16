	CHARACTER (len=7), PARAMETER, PUBLIC :: version='V3.2.12'

!Signification du numéro de version : Va.b.c
!a augmente si les fichiers binaires de reprise *.pms, *.rep, *.dat... changent
!même si certains restent compatibles
!b augmente si l'un des fichiers, données (*.don), réglage, sortie (*.osc)...
!change
!c augmente s'il y a des modifications d'algorithmes (bug, nouvelles
!implantations, ...)

!Journal des aménagements de CESAM2k

!15/11/08 CESAM2k.V3.2.12
! 1 - bsp_gal.f, newspl_gal.f : optimisation, introduction de mg
! 2 - eq_ini_rota4.f, initialise_rota4 : initialisation de U par eq. diff. lin.
! 3 - eq_diff_rota4.f : condition limite externe sur Tau au lieu de U=0
! 4 - coeff_rota4.f, eq_diff_rota4.f : multiplication par nu des coefficients de
!     l'équation de U
! 5 - remplacement des routines conv_cm_reza.f et conv_cgm_reza.f, (corrections
!     communiquées par Reza Samadi)
! 6 - fcmax.f: le cas precision = 'mx' n'était pas envisagé
! 7 - eq_diff_chim.f: mcz_ext --> nuzc_ext
! 8 - evol.f: déplacement des définitions de nuzc_ext et mzc_ext,
!     suppression des "bricolages techniques" pour la base de la rotation
! 9 - mod_atm: extension des tables de Marcs introduction/modifications de
!     trho_4000.f, marcs.f, mod_atm.f (B.Pichon), correction d'une coquille
! 10 - pertm*.f: reformulation des pertes/gains de masse, pertm_waldrom.f
!      correction d'un (petit) bug dans la détermination de la luminosité.
! 11 - newspl_gal.f: appels à bsp1dn limités à l'intervalle [x(1),x(nx)],
!      augmentation du nombre d'abscisses pour intégration de gauss
! 12 - ppcno10.f : Li7 est hors équilibre
! 13 - diffus.f: suppression des "bricolages"
! 14 - eq_diff_rota4.f : optimisation de la condition limite de Phi en surface
! 15 - cesam.f : aménagement des commentaires de la logique de baratine,
!      grille_fixe est imposé avec diffusion du moment angulaire
!      Il est possible de poursuivre un calcul avec des résaux nucléaires
!      différents si les éléments chimiques pris en compte sont identiques Ex
!      Exemple: ppcno9 et ppcno3a9. Amélioration du calcul des pourcentages
!      d'énergie PP, CNO..GRAV, incidences pour les écritures et dans list.f
! 16 - lim_zc.f : adaptation du changement de grille fixe
! 17 - abon_ini.f : introduction des mixtures photosphérique et météoritique
!      de Asplund Grevesse Sauval 05, table crée par Yveline Lebreton
! 18 - calib2k_zams.f, calib2k_pms.f introduction de la possibilité de modifier
!      le Z/X de calibration
! 19 - implantation de diffw_toul.f : coefficients de diffusion simplifiés pour
!      la diffusion du moment cinétique
! 20 - marcs.f, roger.f : écriture de la valeur utilisée de [Fe/H]
! 21 - Aménagements dans la notice
! 22 - module mod_exploit : des2k_abon.f introduction du dessin de grad_mu,
!      création de des2k_dvson.f
! 23 - Suppression des routines ctes*.f, les initialisations sont effectuées
!      directement dans ini_ctes. Les noms des ensembles de constantes
!      i.e. ctes_94 restent inchangés. Création d'ensembles de constantes
!      ctes_31 (de Toulouse), ctes_ba (Basu & Antia) et ctes_gaia.
! 24 - lit_nl.f : possibilité de faire le calcul avec rot_0 "ie.sans rotation"
!      bien que la vitesse angulaire initiale w_rot soit (par erreur) non nulle
! 25 - Introduction de la correction de Jorgen (astro-ph 0811.100v1) coder
!      NOM_OPA=opa_yveline_jorgen
! 26 - write_nl.f : utilisation de grad_ovi/s pour le signe de ovshti/s
! 27 - Atmosphères roger & marcs : utilisation des valeurs limites en cas de
!      sorties de table en teff et/ou logg, adaptations dans tdetau.f,
!      trho.f, trho_4000.f
! 28 - Introduction de tacho_ext/full.f: diffusion sous/sur et sous les ZC 
!      selon les Toulousains, pour Li7 solaire calibration, utilisation avec
!      difft_nut_ext/full 

!20/09/08 CESAM2k.V3.2.11
! 1 - déplacement de coeff_vth de mod_etat --> mod_evol
! 2 - diffw_mps.f prise en compte du statut OPTIONAL de dhv
! 3 - coeff_rota4.f correction sur chi_t
! 4 - coeff_rota.f, coeff_rota4.f, eq_diff_rota4.f suppression des dérivées
!     dhv et dfrl inutilisables
! 5 - lit_nl.f : suppression de l'utilisation de la condition limite JPZ
!     pour la rotation (il conviendrait de passer à la version V3.3...)
! 6 - mise à jour des SCRIPTS_BASH

!09/08/08 CESAM2k.V3.2.10
! 1 - ecrit_ascii.f zbar --> z_bar ligne 366
! 2 - bsp_dis.f n --> nx ligne 81
! 3 - cesam.f psi0 : 0.07 -->0.08 par défaut, 
!    Krot=3,4 : l_demi=.FALSE., ini0=3
! 4 - coeff_rota.f : dfrl(Krot,.. --> dfrl(3,..
! 5 - diffw.f dhv en OPTIONAL, coeff_rota3 suppression de dhv
! 6 - eq_diff_rota4: equation de psi dans les ZC
! 7 - evol : simplification de la gestion des ZC/ZR

!14/07/08 CESAM2k.V3.2.9
! 1 - cesam.f : aménagement des réglages pour la précision 'av'
! proportions d'énergies en % de l'énergie nucléaire totale
! ajustement des réglages avec Krot=3,4
! 2 - evol.f : amélioration de l'approche de la fin de la MS, 
!  introduction de l_demi dans les précisions (cesam.f, mod_donnees.f)
!  amélioration de l'algorithme imposant le mélange radiatif de l'enveloppe
! 3 - lim_zc.f : aménagement de la suppression des limites ZR/ZC trop proches
! 4 - diffw.f : suppression de deff, dv, dh >= 50

!06/05/08 CESAM2k.V3.2.8
!1 - list : attribution de SAVE à fes, bug signalé par B.Pichon
!2 - cesam : introduction de la possibilité de poursuivre avec diffusion
! un modèle initialement sans diffusion et réciproquement (souhait de
! F.Thévenin). Limitation du pas temporel maximal suivant la valeur de la masse
! initiale. Introduction de l_demi. Création du fichier ASCII 4d-2.pms pour
! initialisation PMS de modèles M>10Msol, et ZAMS m=15Msol, le fichier ASCII
! 4d-2.pms et m150.zams sont mis dans le sous directory EXPLOIT
!3 - routines de réactions thermonucléaires : suppression de rot_solid inutile
!4 - eq_diff_chim : réintroduction de Deff avec rotation, voir commentaires
!5 - lim_atm : adjonction d'une loi P(tau) pour Teff > 12000
!6 - mod_donnees : augmentation de pnzc=10. Introduction de l_demi
!7 - evol : introduction de demi, simplification de la logique de suppression
! d'une ZR externe, déplacement de la logique des mvt_dis
!8 - opa.f : la limite pour les opacités compton est ramenée à de 7 à 6.D7
!9 - upddate : suppression de la limitation du pas temporel qui existe dans evol

!01/05/08 CESAM2k.V3.2.7
!1 - etat_opalX/Z : correction de dut=duro*drot+eos(4)*1.e6 au lieu de eos(8)
!    erreur signalée par Yveline, Joséphina, Laurent
!2 - ecrit_ascii correction d'un bug signalé par Yveline. Plantage lors de la
! formation d'un fichier ASCII pour un modèle totalemnt convectif avec ajout
! d'un point au voisinage du centre
!3 - coeff_vth :  adaptation au cas en_masse=.FALSE.

!25/04/08 CESAM2k.V3.2.6
!1 - cesam.f : prise en compte de blabla
!              aménagement des paramètres pour modèles en eulérien
!2 - static_r : remise en fonction
!3 - lim_zc : correction d'un bug affectant la suppression des noyaux radiatifs
!             ou convectifs d'étendue insignifiante

!20/03/08 CESAM2k.V3.2.5
!1 - ascii.f : correction d'un bug signalé par B.Pichon :
! ligne 144 lire evar(2,:)=EXP(evar(2,:)) et non evar(2,:)=10.d0**evar(2,:)
!2 - le fichier journal change de nom, il devient journal.f

!10/02/08 CESAM2k.V3.2.4
!1 - Rotation : rassemblement des routines resout_rota3/4 dans resout_rota.f.
! Création de la routine jacob_rota.f; conditions limites identiques pour les
! formalismes tz97 et mz04; mise en place de l'option lim_jpz pour les deux
! formalismes. Introduction de ord_rot=m_rot+r_qs pour diffusion du moment
! cinétique, ord_rot=m_rot sinon; adaptation de lit_binaire.f
!2 - Suppression de la possibilité d'utiliser les anciens fichiers de données.
!3 - Mise à jour du MAKEFILE du sous-directory EXPLOIT
!4 - Correction d'un commentaire dans marcs.f (cf. B.Pichon).

!01/01/08 CESAM2k.V3.2.3
!1 - cesam.f : Création de la précision 'aj' (ajuste) précision réaliste ('pr')
! avec ajustement du dernier pas temporel pour obtenir  Tc, He core, Xc à
! la fin d'une évolution.
!2 - mod_numerique : création de la routine bsp_gal, permettant la
! détermination des coefficients splines d'interpolation de n fonctions connues
! en tout point. Le but étant d'obtenir une estimation des dérivées.
! Création de la routine de changement de base newspl_gal.
! Ces routines utilisent le formalisme des éléments finis Galerkin. 
!3 - mod_evol : création de la routine tab_vth, tabulation de grandeurs
! thermodynamiques sur le base de la composition chimique. Le but étant le
! calcul des dérivées spatiales de ro, mu.. d'où vaissala. Création de la
! routine coeff_vth permettant le calcul des quantités à tabuler.
! Les tableaux vth et vth_t sont des éléménts de mod_variables.
! Allocation, exploitation de ces tableaux à divers endroits : cesam,
! ecrit_ascii, evol, coeff_rota4.
! 4 - Aménagements pour l'atmosphère. Introduction des lois t(tau) hsra et marcs
! ces dernière sont disponibles sur demande auprès de B.Pichon). Correction
! de bugs dans trho et roger, améliorations de fesh+00.data et fesh+02.data
!(B.Pichon)
! lim_atm : rétablissement de la limite tau_max=20 avec loi T(tau) non
! purement radiative

!01/10/07 CESAM2k.V3.2.2
!1 - Implantation de la routine mu_mol (Module mod_etat) de calcul de diverses
! quantités liées à µ. Introduction dans la routine dgrad pour le critère de
! Ledoux, et implantation dans thermo, lim_zc, for037. Introduction de mu_mol
! dans le MAKEFILE du directory EXPLOIT.
!2 - Amélioration de la définition de phi, incidences dans coeff_rota3/4,
! eq_diff_chim, cesam. et ecrit_ascii sur la fréquence de BV avec new_bv
! Utilisation de new_bv=TRUE pour toutes les précisions sauf 'sa'.
! Dans les fichiers de sortie ASCII, s'il n'est pas tenu compte de la pression
! turbulente, replacement de Ptot/Pgaz=1 [var(23,i)], par grad_mu=dln mu/dln P.
!3 - Le critère de Ledoux peut être utilisé avec et sans convection.
!4 - Adjonction de g, msol, rsol et lsol dans les glob des fichiers ASCII de
! sortie.
!5 - Dans lim_zc : réduction du nombre d'allocations, déplacement de la
! déallocation de new toujours alloué (ligne 1157,remarque de B.Pichon).
! Déplacement du calcul et des allocations de v_son et r_son 
! (remarque de B.Pichon).
!6 - Implémentation de la semi-convection : routine difft_smc du type difft.
! Introduction dans le MAKEFILE du directory EXPLOIT.
!7 - Dans dnunl suppression des allocations (remarque de B.Pichon).
!8 - Dans cesam, update, (reglages) introduction de la variable de "précision"
! dlntc limitant la variation de la température centrale sur un pas temporel,
! Correction d'une erreur pour le calcul BV dans la cas new_bv=.TRUE..
!9 - Dans lit_nl : suppression de la possibilité de l'utilisation des anciens
! fichiers de données (suggestion de B. Pichon).
!10 - Sous-directory EXPLOIT création du programme de dessin des2k_grad,
! des2k_diff_spl, amélioration de f037_2k, et de divers programmes de dessin.
!11 - Mise à jour de la NOTICE.
!12 - Sous-directory SCRIPTS, création de scripts en BASH.

!14/07/07 CESAM2k.V3.2.1
!Dans cesam :
!1- Avec diffusion du moment cinétique (Krot=3, 4) afin de laisser deff jouer
! pleinement son rôle, re_nu est fixé à 0
!2- Introduction du facteur fmin_abon, défini par la précison ou les réglages,
! fmin_abon permet de régler ab_min dans le but d'une amélioration
! du pas temporel lorsque H brûle en couche. Précédemment
! fmin_abon=0.01 était fixé.
!3- Dans ecrit_ascii, ajout de points dans les couches centrales, 
! localisation suivant les paramètres q0 et l0
!4- Dans mod_evol, evol : introduction de coll_rot et colpnt_rot pour la
! détermination des points de collocation pour la rotation. Un décentrage des
! points au milieu des l'intervalles permet de stabiliser la solution avec
! m_rot=2, 4,..pair
!5- Dans EXPLOIT : introduction de coll_rot et colpnt dans le MAKEFILE
!6- Dans evol : amélioration de l'algorithme de l'adaptation des la ZC externe.
!7- Dans coeff_rota3/4 phi_e ne peut être nul, + schéma implicite

!01/07/07 CESAM2k.V3.2.0
!Dans le fichier reglages : ajout de l0 points à disposer autour des
! discontinuités dans les fichiers de sortie ASCII, ajout de la variable
! new_bv calcul de la fréquence BV utilisant Phi (d ln ro / d ln mu),
! avec, faute de mieux, Phi=1
!Dans ecrit_ascii : adaptation pour ajout des points
!Dans bsp_dis et noeud_dis, adaptation au formalisme de cesam2k 
! du traitement des discontinuités dans les routines
!Dans evol : s'il y a retrait du coeur convectif
! 1- sans diffusion, interpolation linéaire de la comp.chim. dans la
! la zone de retrait, 
! 2- avec et sans diffusion suppression de la discontinuité au niveau du
! raccord ZC/ZR.
! Mise en réserve du traitement de l'augmentation de l'abscisse inférieure
! des ZC externes.

!Dans cesam : m_ch=2 (interpolation linéaire) de la composition chimique
! pour toutes les précisions
!Dans mod_numerique implantation de max_local : détermination des maxima locaux
! pour normalisation dans les dessins
!Dans EXPLOIT, modification des programmes des2k_vaiss et des2k_abon, mise en
! place de max_local dans le MAKEFILE

!14/06/07 CESAM2k.V3.1.3
!Ajout de la routine ctes_94m, identique à ctes_94 avec des valeurs entières
!pour les masses atomiques.
!Dans la routine ecrit_ascii, ajout éventuel d'un point à une distance q0 du
!centre. Ajout de q0 dans le fichier reglage.
!Dans le sous-directory EXPLOIT :
! Installation du zoom dans le programme des2k_abon.
! Modification du zoom dans le programme des2k_vaiss.

!25/05/07 CESAM2k.V3.0.2
!Dans eq_diff_rota3/4, resout_rota3/4 :
! aménagement des conditions limites
!Dans sous-directory EXPLOIT :
! update des jacobiens pour les calibrations en tenant compte de la diffusion du
! moment cinétique, programmes calib2k_pms, calib2k_zams
!Dans le sous-directory NOTICE :
! correction de quelques points, des conditions limites pour la rotation

!03/05/07 CESAM2k.V3.0.1
!Dans cesam.f :
! correction d'un bug concernant les proportions pp, CNO, 3a etc..
! m_rot=2 pour la précision 'np'
!Dans coeff_rota4 : élimination du cas nui <= 0
!Dans coeff_rota3 : introduction de phi/delta*gad_mu dans dgrad (remarque de Andy)
!Dans eq_diff_rota3 et eq_diff_rota4 :
! correction d'un bug concernant la dérivée de bs(2)
! calcul direct de Lambda, reprenant (idée de Phi)
!Mise en place de la variable logique baratine permettant de
! détourner une grande partie des commentaires 'on line' vers les fichiers
! mon_modele_static, *_atmos, *_evol pour, respectivement, l'équilibre
! quasi-statique, l'atmosphère, l'évolution de la composition chimique et
! la rotation (idée de B.Pichon).
!Dans lit_nl introduction de la possibilité de lire un fichier .don
! ancien et un fichier .don sans diffusion ni rotation (idée de B.Pichon).
!Création du fichier pgplot_factice.f permettant d'éviter l'utilisation du
! PGPLOT pour les dessins "on line"

!03/05/07 CESAM2k.V3.0.0
!Pour la rotation, remplacement de la méthode des éléments finis par celle
!de collocation
!Pour des initialisations, les fichiers binaires *.hom, *.pms sont utilisables

!03/02/07 CESAM2k.V2.5.2
!Pour la rotation : mise en place de formalismes identiques dans ZR et ZC
!avec, dans ZC, des coefficients de diffusion >> 1

!04/02/07 CESAM2k.V2.5.1
!Modification du fichier reglages
!Dans evol.f :
! Suppression de demi
!Dans cesam.f :
! Avec Krot=3,4 on ignore les discontinuités de composition
! chimique et on les lisse mvt_dis=.FALSE., lisse = .TRUE.
! pour la diffusion du moment cinétique utilisation de la formulation approchée
! de mu(Krot=3,4), mu_saha=.FALSE. (sauf avec precision='rg')
! Pour précision réaliste & corot (pr & co),
! fonction d'espacement limitée à Ln P, ajustement de psi0 pour avoir un nb de
! couches de l'ordre de 1000 pour une 1.5Msol sur la ZAMS
! addition de ajuste et lisse dans le fichier reglages
!Dans base_chim.f :
! Utilisation d'une base continue non dérivable
!Dans base_rota.f :
! Utilisation d'une base continue non dérivable
!Dans print_ctes.f :
! permutation de Li6 et Li7 dans la liste d'écriture
!Dans coeff_rota3.f, coeff_rota4.f, eq_diff_poisson.f :
! allocation du tableau ion, (cf. B.Pichon)

!19/01/07 CESAM2k.V2.4.6
!Dans coeff_rota4
! Addition de pi dans cte5_0
! Changement de signe de cte_0(4) Bug signalé par Yveline
! Corrections des dln_mu Bug signalé par ANDY
! Addition de mu pour les tracés :
!  Augmentation à 30 de nb_coeff dans mod_evol
!  Modification de des2k_rot dans EXPLOIT
!Dans resout_rota4, ligne 345 2 au lieu de iv

!11/12/06 CESAM2k.V2.4.5
!Dans coeff_rota3 et 4
!Introduction d'une approximation de dln ro / dln mu
!Utilisation de l'approximation numérique de d ln P/d nu
!Dans cesam.f, énergies pp, cno, 3a+C+O, grav en Lsol
!Dans evol permutation de l'ordre des intégrations de diffusion
!composition chimique puis moment cinétique
!Dans coeff_rota4, utilisation de chi_T
!Dans resout_rota4 suppression du controle des corrections NR

!06/12/06 CESAM2k.V2.4.4
!Aménagement des tests de dérivation dans coeff_rota4, eqdiff_rota4
!Correction d'Yveline : dans cesam.f, utilisation de -TdS au lieu de ABS(TdS)
!pour l'estimation de l'énergie graviphique.

!25/11/06 CESAM2k.V2.4.3
!Amélioration de la logique du dtmax dans cesam.f
!Suppression des discontinuités dans la base de la rotation
!les variables deviennent continues non dérivable ==> suppression du
!rétablissement des continuités dans resout_rota3/4.
!Augmentation du nombre de variables dans le fichier ASCII *_coeff-rota.dat
!Aménagement du programme des2k_rot du sous directory EXPLOIT pour le tracé des
!des variables de la rotation
!Limitation à 50My du pas temporel s'il y a diffusion du moment cinétique
!et limitation à 10% la variation du pas temporel
!Augmentation/diminution du nombre de couches limité à 5% (au lieu de 10%)

!27/10/06 CESAM2k.V2.4.2
!rectification lnt923" --> lnt932  réaction 5 de NACRE dans taux_nuc.f
!rectification des commentaires concernant mu_e, var(15) est mu_e
!Mise à jour de la notice et de l'aide mémoire

!12/09/06 CESAM2k.V2.4.1
!Introduction de la variable logique nc_max permettant d'imposer le nombre
!maximum de couches n_max pour le calcul du dernier modèle : nc_max=n_max < 0
!Dans cesam.f, introduction de la variable coox et du fichier *.coox
!combustion de l'oxygène
!Dans opa_houdek9 appel à opa_opal2_co en cas de sortie de table (pour
!l'atmosphère opa_opal2_co, utilise opa_yveline)

!22/08/06 CESAM2k.V2.4.0
!Introduction des réactions nucléaires de la combustion du carbone
!C12(C12,..).. C12(O16,..) et de l'oxygène O16(O16,..).
!Intervalle de tabulation des taux des réactions fixé à 1MK
!Création de la précision 'av' avec des aménagement permettant d'atteindre
!les stades avancés
!Création des routines ppcno3aco  ppcno3acos de la combustion de H à O
!Introduction du vecteur iter_qs permettant d'adapter les variables contrôlées
!pour la résolution de l'équilibre quasi-statique
!Création de la routine opa_compton et utilisation dès que T > 0.007T9
!Ecriture de Teff dans le dessin
!Calcul des poids statistiques dans saha
!Dans des_m abondances centrales avec 'av'

!05/07/06 CESAM2k.V2.3.3
!Introduction des réactions nucléaires de NACRE pour le 3alpha et C12(a,g)O16
!Début de la combustion du carbone mise à 6d8K, tabulation pour ppcno3ac10
!repousée à T9
!Création de la précision 'av' : long runs pour les longues évolutions
!Réaménagement du module mod_nuc
!Mise en place de SAVE dans etat_opal* (remarque de JP Marques)

!28/06/06 CESAM2k.V2.3.2
!Modification de la localisation des limites ZR/ZC dans lim_zc
!Pour le calcul des coefficients de diffusion du moment cinétique, utilisation
!de Omega et U au temps t

!13/06/06 CESAM2k.V2.3.1
!Adaptation de l'expression approchée de \phi
!Possibilités de dessin on line et off line des coefficients de la diffusion
!du moment cinétique

!24/05/06 CESAM2k.V2.3.0
!Introduction du formalisme de diffusion du moment cinétique
!selon Mathis & Zahn 2004 occasionnant diverses adaptations du fichier de
!données et de réglages.
!Création d'une routine de dessin des variables de la diffusion du moment
!cinétique
!Utilisation de l'ancienne formule de la fréquence de Brunt Vaissala
!Correction d'un bug dans les routines PPCNO12, PPCNO12Be, PPCNO12BeBFe,
!PPCNO12Li
!Accélération de la convergence de la diffusion du moment cinétique

!24/04/06 CESAM2k.V2.2.0
!Restructuration permettant différents calculs de la vitesse angulaire
!Modification du fichier de données
!Introduction du formalisme de Matis & Zahn 2004 (début)
!Introduction de la conservation locale du moment cinétique

!11/04/06 CESAM2k.V2.1.0
!Vérification du Jacobien de ppcno3a12Ne (bug)
!Redéfinition de scale dans rkimps
!Suppression, dans cesam, de Kipp=.TRUE. pour les modèles après la séquence
!principale
!Suppression, dans cesam, de n_max de la NAMELIST nl_rlg

!20/03/06 CESAM2k.V2.0.8
!Aménagements dans etat_opal, ZFSinterppeos, opal_ascii_bin, calib2k
!Suppression des SUM dans opa_opal2 (bug)
!Compléments de formules de réactions nucléaires NACRE
!Aménagements mineurs dans des_m, des_r, resout, z14xcotrin21

!20/03/06 CESAM2k.V2.0.7
!Création du programme d'exploitation des2k_opa.f
!Création de la routine générique coeff_rota.f et des routines coeff_rota_saha.f,
!coeff_rota_ioni.f, coeff_rota_z16.f
!Dans ces routines mise à 0 de chi_T mal calculé avec les données dont on dispose
!Création de la routine difft_sun.f
!Addition de la viscosité cinématique au coefficient Dv
!Création de la routine coeff_rota_ioni.f
!Introduction de mini, valeur minimale de Dv, dans le module mod_evol
!Dans lit_nl, avec diffusion du moment cinétique, on impose D_turb >= 50
!Introduction du nom de la routine de calcul des coefficients de rotation dans
!le type de précision et dans le fichier réglages
!Introduction dans mod_donnees de la variable logique ecrit_rot conditionnant
!l'écriture du fichier mon_modele_coeff_rota.dat pour dessin des coefficients
!de rotation, addition de ecrit_rot dans le fichier reglages
!Création de la routine difft_gab.f
!Mise à jour des paramètres de précision et de l'aide mémoire
!Aménagements du test de sécurité MODIF_CHIM dans abon_ini.f vent.f,
!planetoides.f

!02/02/06 CESAM2k.V2.0.6
!Suppression de w_form du fichier de données, mis dans le fichier reglage
!Création de la routine générique coeff_rota.f appelant l'une des deux routines
!coeff_rota_saha.f et coeff_rota_z16.f
!Mise en place des chutes de planétoïdes : modification de la composition
!chimique de la ZC externe et apport/retrait de moment cinétique, adaptation du
!fichier planet

!01/02/06 CESAM2k.V2.0.5
!Limitation de la source du vent à la ZC externe, suppression du
!paramètre p_vent du fichier vent, simplification du traitement du vent,
!aménagement du programme fichier_vent.f et des fichiers exemple.vent et vent
!Mise en place des chutes de planétoïdes, création de la routine planétoïdes et
!des fichiers exemple.planet et planet
!Ajonction de la masse terrestre dans les fichiers ctes85 et ctes94
!Ajonction de W_FORM dans le fichier "reglages"
!En cours :
!Création, suppression, aménagement et mise au point des divers routines et
!programmes concernés par la diffusion du moment cinétique, principalement :
!cesam, ecrit_rota, diffus,
!resout_rota, resout_chim, coeff_rota, eq_diff_rota, eq_diff_chim, des2k_rot

!15/12/05 CESAM2k.V2.0.4
!Introduction de l'argument optionnel duale de la routine newspl
!Construction et aménagements de diverses routines pour la diffusion du moment
!cinétique

!20/10/05 CESAM2k.V2.0.3
!Aménagements mineurs dans inter
!Correction du calcul de d ln l / d ln m, et d ln ro / d m^2/3
!Possibilité de calcul direct des coefficients de diffusion du moment cinétique
!Introduction de tab_coeff_rota
!Introduction dans resout de la variable et fonction logique cmax et fcmax pour
!l'utilisation du nombre maximum de couches avant de sortir
!Suppression du dessin de ro sans diffusion (discontinuités)
!Dans evol, sans diffusion,léger lissage par contour de la composition chimique
!pour lisser le retrait des ZC
!Ajout de lim et de model_num à la fin des fichiesr binaires

!13/10/05 CESAM2k.V2.0.2
!Calcul direct de Deff dans eq_diff_chim
!Permutation de l'ordre diffusion du moment cinétique <==> diffusion des
!éléments chimiques
!Création du programme de dessin des2k_coeff_rota
!Ecriture du numéro du modèle dans des_m et des_r
!Dans les fichiers de sortie ASCII, ajout de 20 points de grille de part et
!d'autre des limites ZR/ZC pour affiner le profil de la fréq.BV
!En abscence de diffusion amélioration de la formulation de la fréq.BV

!05/10/05 CESAM2k.V2.0.1
!Amélioration du choix de no_croiss dans noein
!Amélioration d'écritures et introduction de no_croiss dans linf
!SAVE pour les quantités ***0 des conditions limites de static_m

!01/10/05 CESAM2k.V2.0.0
!SAVE dans opa_yveline, etat_opalX, etat_opalZ
!Mise en place des numéros des modèles, sorties de tous les modèles en ASCII
!et .rep avec leur numéro, conservation du numéro dans les fichiers binaires
!Choix de grad_ad ou grad_rad dans les zones overshootées

!20/09/05 CESAM2k.V1.1.15
!Adjonction de v dans l'expression de teta dans les tests de dérivation de
!static_m et static_r
!Création de la routine ppcno3a13Ne22

!01/09:05 CESAM2k.V1.1.14
!Correction de dgravr dans thermo et thermo_atm
!Facteur 2/3 sur l'accélération centrifuge dans coll_atm et eq_atm
!Aménagements des équations relatives à la diffusion du moment cinétique
!Création des programmes de dessin des2k_dhve, des2k_rot, des2k_bin

!30/08/05 CESAM2k.V1.1.13
!Suppression de commentaires dans les modules mis dans la notice
!Suppression de nom_elem en dp des définition de mod_nuc
!Suppression de la variable pmw du module mod_donnees
!Suppression des tableaux xlim_rot, xcin et xcint du module mod_evol

!03/08/05 CESAM2k.V1.1.12
!Inversion de la chronologie du journal
!Corrections de bugs signalés par A.Moya dans coeff_rota :
!C12=1 et signes - pour C16,17,18
!Rétablissement de d2U/dnu2=0, et équation de diffusion de Omega dans ZC
!Déplacement de l'allocation de frot, dfrot dans lim_zc
!Annulation de C15*, C8 et C9 dans coeff_rota

!27/06/05 CESAM2k.V1.1.11
! coeff_rota, utilisation de rho ie. sans passer par l'équation d'état,
! pour cohérence avec dln ro
! des_m, augmentation du nb. de chiffres significatifs pour les abondances max
! diffm_mp changement de signe de l'accélération centrifuge
! diffm_mp & diffm_br coefficient de l'accélération centrifuge
! aménagements dans resout_rota, eqdiff_rota, diffus, coef_rota
! création du programme de dessin des2k_dhve du sous-directory EXPLOIT

!16/06/05 CESAM2k.V1.1.10
! Suppression du fichier *.atm pour initialiser ZAMS ou PMS
! Tracé de ro dans des_m et des_r
! Abondances des éléments au centre dans list
! Augmentation de m_rot --> 4 dans cesam.f pour tous les réglages
! Après la ZAMS on impose l'approximation de Kippenhahan

!13/06/05 CESAM2k.V1.1.9
! Correction C12(a,g)O16 dans ppcno3ac10
! Elargissement des dessins du HR jusqu'à nb nmax modèles
! etat_opal SAVE pour la variable iri
! Lxchim(nchim) dans etat_opal et etat_ceff
! lit_nl, nb_max_models pour lit_nl_2korg
! Réduction du pas temporel à la fin de la ZAMS
! Augmentation du nombre de couches après la TAMS, He et C burning
! Définition de z_table=z0 dans opa_opal2
! diffw routine PUBLIC de mod_evol
! Création de lit_binaire dans mod_exploit
! Création du programme des2k_dhve dans EXPLOIT

!03/06/05 CESAM2k.V1.1.8
! Rectification d'une virgule dans mod_exploit
! Suppression de la référence à compg(ihe4,1) dans list
! Save de cte1 dans colatm
! Mise du numéro de version dans journal

!31/05/05 CESAM2k.V1.1.7
! Mise en service du programme des2k_rot du sous-directory EXPLOIT
! SAVE et allocation des tables de données dans opa_yveline

!26/05/05 CESAM2k.V1.1.6
! Restriction de l'utilisation de ln ro = bp(7,:) au cas avec diffusion
! et ord_qs  2 à cause de la discontinuité de ro aux limites ZR/ZC
! Introduction du numéro de version par un include dans mod_donnees

!20/05/05 CESAM2k.V1.1.5
! Corrections dans lim_gong1, lim_tau1, add_ascii, opa_yveline lisse,
! pour réimplantation des modèles de GONG

!14/05/05 CESAM2k.V1.1.4
!Sous directory SOURCE:
! Correction du bug dans la formule de normalisation dans evol
! lire	chim(1:nchim,i)=chim(1:nchim,i)/norm
! et non  chim(1:nchim,:)=chim(1:nchim,:)/norm
!Sous directory EXPLOIT:
! Création des programmes de dessin des2k_abonts et des2k_abontc
! Suppression du programme des2k_abont
! Dessin de X, Y, Z dans des2k_abon

!05/05/05 CESAM2k.V1.1.3
! Correction d'un bug lié au calcul de Teff si n_atm=1 ie. dans les
! cas GONG1 et GONG2
! Calcul ("exact") de la féquence de Brunt-Vaissala en utilisant ln ro
! Suppression de contour

!02/05/05 CESAM2k.V1.1.2
! Amélioration des algorithmes gérant les arrêts sur t_stop et x_stop
! Implantation de l'arrêt sur He_core

!27/04/05 CESAM2k.V1.1.1
!  Correction de bugs engendrés par ihe4=-100 et lvent=.TRUE.avec PP1

!23/04/05 CESAM2k.V1.1.0
! Introduction dans resout de l'arrêt sur x_stop
! Introduction dans resout de l'arrêt sur t_stop

