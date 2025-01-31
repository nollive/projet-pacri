// Permet de lancer des simulations en parallele a partir d une grille de parametre qui va etre genere
// TUTO
// pour faire tourner la pipeline nextflow, il faut se placer dans un dossier et avoir dedans ce script, le script nextflow.config, et un dossier bin avec les codes R
// il faut également le fichier code_model.R et renseigner son path dans ce fichier
// on lance la ligne de commande suivante : nextflow run main.nf --simu_count 2
// parametre simu_count : nombre de simulations effectuees avec un jeu de parametre

// pour modifier les valeurs des parametres de la grille, il faut modifier le fichier generation_param.R
// pour modifier les parametres à faire varier, il faut modifier les fichiers generation_param.R, simu_nextflow.R et stats_simu.R
// pour modifier les indicateurs des simulations, il faut modifier le fichier simu_nextflow.R et stats_simu.R

// pour utiliser sur le cluster de Pasteur, il faut d abord charger les modules avec : module load apptainer, module load graalvm et module load nextflow


nextflow.enable.dsl=2

// parametres
params.simu_count = 2
params.folder_result = 'results'
params.folder_simu = 'simu'

// chemin entier vers le code R du modele
// sur le cluster
// params.rlib= "/pasteur/zeus/projets/p02/Pacri/maylayan/Hospitals_modes_transmission/github_codes/cpp/dev-sensibility-analysis.cpp"
// params.rlib2= "/pasteur/zeus/projets/p02/Pacri/maylayan/Hospitals_modes_transmission/github_codes/R/nodscov2/helper-functions-simulations.R"

// en local
params.rlib="/Users/maylayan/Documents/Projets/Hospitals_modes_transmission/github_codes/cpp/dev-sensibility-analysis.cpp"
params.rlib2="/Users/maylayan/Documents/Projets/Hospitals_modes_transmission/github_codes/R/nodscov2/helper-functions-simulations.R"

// chemin vers le dossier ou se trouve les fichiers .RData input
// params.rinput = "/pasteur/zeus/projets/p02/Pacri/maylayan/Hospitals_modes_transmission/github_codes/data/data-synthetic-graphs/loc/"
params.rinput = "/Users/maylayan/Documents/Projets/Hospitals_modes_transmission/github_codes/data/data-synthetic-graphs/loc/"




process generation_param {
	// permet de generer la grille de parametre a partir du nombre de simulation que l on veut effectuer pour chaque jeu de parametre
    publishDir "${params.folder_result}", mode: 'copy' //sous dossier de Nextflow ou seront copies les output 

    input:
	  val simu_count //valeur du nombre de simulation a realiser pour chaque jeu de parametres

    output:
	  path 'param_grid.txt' //fichier .txt avec une ligne par set de parametres, chaque colonne represente un parametre

    script:
    """
	  Rscript ${workflow.projectDir}/bin/generation_param.R ${simu_count}
    """
}

process simulation {
	// realise une simulation a partir d un jeu de parametre et calcul des indicateurs interessants
	publishDir "${params.folder_result}/${params.folder_simu}", mode: 'copy'
  
  // queue { 
  //  if (task.attempt < 3) {
  //    'common,dedicated'
  //  }  else if (task.attempt < 4) {
  //    'common'
  //  } else {
  //    'long'
  //  }
  // }
  
  clusterOptions { 
    if (task.attempt < 3) {
      '--qos=fast'
    }  else if (task.attempt < 4) {
      '--qos=normal'
    } else {
      '--qos=long'
    }
  }

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  
	input:
	path rlib //chemin d acces au code R helper-functions-cpp.R
	path rlib2
	path rinput //chemin du dossier avec les differents input.Rdata
  val line //une ligne du fichier param_grid.txt, qui correspond a un vecteur avec un set de parametres
	
	output:
  path 'out*.rda' //fichier RData qui sort du modele 
	path 'summary_stat*.csv', emit: resu //fichier .txt avec une seule ligne avec tous les indicateurs de la simulation

  script:
  """
	Rscript ${workflow.projectDir}/bin/simu_nextflow.R $rlib $rinput ${line}
  """
}

workflow {
	// genere la grille de set de parametre 
  param_grid = generation_param(params.simu_count) 
	
	// permet de split chaque ligne de param_grid 
	lines = param_grid
		.splitText(keepHeader: true)

	// fait tourner la simulation sur chaque set de parametre
	rlib=file(params.rlib)
	rlib2=file(params.rlib2)
	rinput=file(params.rinput)
	simulation=simulation(rlib, rlib2, rinput, lines)
	
	// collecter tous les indicateurs de chaque simu dans un seul fichier que l on enregistre dans le dossier de resultat
   simulation.resu
		.collectFile(name: 'resu_simu_all.txt', newLine: true, keepHeader: true, skip: 1)
		.subscribe { file -> file.copyTo("${params.folder_result}/resu_simu_all.txt") }
		
	// utiliser tous les indicateurs de chaque simu pour faire des stats dessus
	//all_simu = simulation.resu
		//.collect()
	//stat_simu(params.simu_count, param_grid, all_simu )
}