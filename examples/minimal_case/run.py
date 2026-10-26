"""Run the minimal case without an external SWAT project."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from swat_dps.optimize import main
if __name__=='__main__':
    defaults=['--config',str(Path(__file__).with_name('config.json')),'--population','20','--generations','5','--workers','2',
              '--output',str(Path('results')/'minimal_case')]
    main(defaults+sys.argv[1:])
