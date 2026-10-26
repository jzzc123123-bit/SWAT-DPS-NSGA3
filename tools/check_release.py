"""Check current source hashes, documentation links and release data boundaries."""
from pathlib import Path
import hashlib,json,re,subprocess

def main():
    root=Path(__file__).resolve().parents[1]
    manifest=json.loads((root/'docs/swat_source_manifest.json').read_text(encoding='utf-8'))
    for name,digest in manifest['sources'].items():
        p=root/'swat_source/dps_b8k3'/name
        if hashlib.sha256(p.read_bytes()).hexdigest()!=digest:raise ValueError('Source snapshot hash differs: '+name)
    for p in [root/'README.md',*(root/'docs').glob('*.md'),root/'examples/minimal_case/README.md']:
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if re.match(r'https?://',target):continue
            target=target.split('#',1)[0]
            if target and not (p.parent/target).exists():raise ValueError('Broken documentation link: '+target)
    files=subprocess.check_output(['git','ls-files','--cached'],cwd=root,text=True).splitlines()
    allowed_roots={'.gitignore','.gitattributes','README.md','requirements.txt','swat_nsga3.py'}
    allowed_docs={'docs/METHOD.md','docs/SWAT_INTEGRATION.md','docs/swat_source_manifest.json'}
    allowed_dirs=('swat_dps/','swat_source/dps_b8k3/','examples/','tools/','tests/','.github/')
    unexpected=[p for p in files if p not in allowed_roots and p not in allowed_docs and not p.startswith(allowed_dirs)]
    if unexpected:raise ValueError('File outside the source/example release: '+unexpected[0])
    if any(p.startswith('RawData/') for p in files):raise ValueError('Real-project RawData remains in the publication index')
    if any(p.startswith(('audit_local/','results/','build/')) for p in files):raise ValueError('Private/generated run files staged')
    for directory in ['swat_dps','examples','tools','tests']:
        for p in (root/directory).rglob('*'):
            if not p.is_file() or '__pycache__' in p.parts:continue
            if p.suffix.lower() not in ['.py','.json','.md','.f90']:raise ValueError('Unexpected release data/binary: '+p.name)
            if re.search(r'(?m)(?:gh[pousr]_[A-Za-z0-9]{20,}|^-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----$)',p.read_text(encoding='utf-8')):
                raise ValueError('Possible credential in '+p.name)
    print(json.dumps({'passed':True,'current_swat_sources':len(manifest['sources']),
                      'research_inputs_in_current_tree':False,'generated_or_private_runs_staged':False}))

if __name__=='__main__':main()
