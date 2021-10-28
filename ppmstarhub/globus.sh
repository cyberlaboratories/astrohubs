echo Function globus_help provides help for Globus CLI

globus_help(){
    echo ""
    echo Quickstart: https://docs.globus.org/cli/quickstart
    echo ""
    echo Start with `globus login --no-local-server` which will guide you through the authentication process.
    echo ""
    echo Some example commands:
    echo globus endpoint search ranch
    echo globus ls d98c7f06-6d04-11e5-ba46-22000b92c6ec:/stornext/ranch_01/ranch/projects/AST20006/PPMstar
    echo or, using the variables, see list below:
    echo globus ls $niagara_scratch
    echo globus mkdir $astrohub__dataA/testdir
    echo globus endpoint search Frontera
    echo globus transfer -r --dry-run $niagara_scratch_fherwig/M111/png $fherwig__laptop_Globus 
    echo globus task list
    echo globus task show 716769e4-1dfb-11eb-81b7-0e2f230cc907
    echo ""
    echo Globus endpoint variables with some often used directories
    echo we will use double underscore for # in endpoint names and single
    echo underscore for forward slash in the directory part
    echo ""
    echo "Here are pre-defined endppoint variables:"
    
    echo computecanada__hermes_scratch_ASDR="b4c53282-6d04-11e5-ba46-22000b92c6ec:scratch/ASDR"
    echo frontera_scratch3_mao="142d715e-8939-11e9-b807-0a37f382de32:/~/scratch3_mao"
    echo frontera_scratch3_paul="142d715e-8939-11e9-b807-0a37f382de32:/~/scratch3_pwoodwar"
    echo ranch_ast20006="d98c7f06-6d04-11e5-ba46-22000b92c6ec:/stornext/ranch_01/ranch/projects/AST20006"
    echo ranch_ppmstar="d98c7f06-6d04-11e5-ba46-22000b92c6ec:/stornext/ranch_01/ranch/projects/AST20006/PPMstar"
    echo niagara_scratch="77506016-4a51-11e8-8f88-0a6d4e044368:/scratch/f/fherwig"
    echo niagara_scratch_fherwig="77506016-4a51-11e8-8f88-0a6d4e044368:/scratch/f/fherwig/fherwig"
    echo fherwig__laptop_Globus="de463cef-6d04-11e5-ba46-22000b92c6ec:Globus"
    echo astrohub__nugrid="13897778-9894-11ea-b3c4-0ae144191ee3"
    echo astrohub__dataA="55dbc7da-80ec-11ea-97a5-0e56c063f437"
    echo astrohub="4af3403c-7f43-11ea-97a5-0e56c063f437"
    echo astrohub__stellarhydro="2c4a17ea-c58f-11ea-8f25-0a21f750d19b"
    echo astrohub__TMT_jet="79a64fc6-80ee-11ea-aff4-0201714f6eab"
    echo astrohub__scratch14_wendi="82f29ba4-9898-11ea-bf91-0e6cccbb0103"
    echo cedar__scratch_fherwig="c99fd40c-5545-11e7-beb6-22000b9a448b:scratch"
    echo hermes__scratch="b4c53282-6d04-11e5-ba46-22000b92c6ec:scratch"
    echo hermes__scratch_ASDR="b4c53282-6d04-11e5-ba46-22000b92c6ec:scratch/ASDR"
    echo niagara__hpss='36002f06-44b1-11e8-8e11-0a6d4e044368:/archive/f/fherwig'
    echo niagara__hpss_hydroArchive='36002f06-44b1-11e8-8e11-0a6d4e044368:/archive/f/fherwig/hydroArchive'
    echo niagara__hpss_otherArchive='36002f06-44b1-11e8-8e11-0a6d4e044368:/archive/f/fherwig/otherArchive'
}
computecanada__hermes_scratch_ASDR="b4c53282-6d04-11e5-ba46-22000b92c6ec:scratch/ASDR"
frontera_scratch3_mao="142d715e-8939-11e9-b807-0a37f382de32:/~/scratch3_mao"
frontera_scratch3_paul="142d715e-8939-11e9-b807-0a37f382de32:/~/scratch3_pwoodwar"
ranch_ast20006="d98c7f06-6d04-11e5-ba46-22000b92c6ec:/stornext/ranch_01/ranch/projects/AST20006"
ranch_ppmstar="d98c7f06-6d04-11e5-ba46-22000b92c6ec:/stornext/ranch_01/ranch/projects/AST20006/PPMstar"
niagara_scratch="77506016-4a51-11e8-8f88-0a6d4e044368:/scratch/f/fherwig"
niagara_scratch_fherwig="77506016-4a51-11e8-8f88-0a6d4e044368:/scratch/f/fherwig/fherwig"
fherwig__laptop_Globus="de463cef-6d04-11e5-ba46-22000b92c6ec:Globus"
astrohub__nugrid="13897778-9894-11ea-b3c4-0ae144191ee3"
astrohub__dataA="55dbc7da-80ec-11ea-97a5-0e56c063f437"
astrohub="4af3403c-7f43-11ea-97a5-0e56c063f437"
astrohub__stellarhydro="2c4a17ea-c58f-11ea-8f25-0a21f750d19b"
astrohub__TMT_jet="79a64fc6-80ee-11ea-aff4-0201714f6eab"
astrohub__scratch14_wendi="82f29ba4-9898-11ea-bf91-0e6cccbb0103"
cedar__scratch_fherwig="c99fd40c-5545-11e7-beb6-22000b9a448b:scratch"
hermes__scratch="b4c53282-6d04-11e5-ba46-22000b92c6ec:scratch"
hermes__scratch_ASDR="b4c53282-6d04-11e5-ba46-22000b92c6ec:scratch/ASDR"
niagara__hpss='36002f06-44b1-11e8-8e11-0a6d4e044368:/archive/f/fherwig'
niagara__hpss_hydroArchive='36002f06-44b1-11e8-8e11-0a6d4e044368:/archive/f/fherwig/hydroArchive'
niagara__hpss_otherArchive='36002f06-44b1-11e8-8e11-0a6d4e044368:/archive/f/fherwig/otherArchive'
