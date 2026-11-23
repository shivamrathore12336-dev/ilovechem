#!/usr/bin/env python3
"""
Agentic Awesome Skills - Fast On-Demand Skill Finder & Loader for Antigravity
Usage:
    python find_skill.py search <keyword>
    python find_skill.py install <skill-id>
    python find_skill.py info <skill-id>
"""

import sys
import json
import shutil
from pathlib import Path

REPO_ROOT = Path(r"C:\Users\shiva\.gemini\antigravity\scratch\agentic-awesome-skills")
INDEX_FILE = REPO_ROOT / "skills_index.json"
SKILLS_DIR = REPO_ROOT / "skills"
GLOBAL_SKILLS = Path(r"C:\Users\shiva\.gemini\config\skills")
WS_SKILLS = Path(r"C:\Users\shiva\.gemini\antigravity\scratch\chemsite\.agents\skills")

def load_index():
    if not INDEX_FILE.exists():
        print(f"Index file not found at {INDEX_FILE}")
        sys.exit(1)
    with open(INDEX_FILE, "r", encoding="utf-8") as f:
        return json.load(f)

def search_skills(query):
    data = load_index()
    q = query.lower()
    matches = []
    for item in data:
        sid = item.get("id", "").lower()
        desc = item.get("description", "").lower()
        cat = item.get("category", "").lower()
        if q in sid or q in desc or q in cat:
            matches.append(item)
    
    print(f"\nFound {len(matches)} skill(s) matching '{query}':\n")
    for m in matches[:15]:
        print(f"  * {m['id']} [{m.get('category', 'unknown')}]: {m.get('description', '')[:95]}...")
    if len(matches) > 15:
        print(f"\n  ... and {len(matches) - 15} more. Refine your query or install using: python find_skill.py install <id>")

def install_skill(skill_id):
    src = SKILLS_DIR / skill_id
    if not src.exists():
        print(f"Error: Skill '{skill_id}' not found in {SKILLS_DIR}")
        return
    
    GLOBAL_SKILLS.mkdir(parents=True, exist_ok=True)
    WS_SKILLS.mkdir(parents=True, exist_ok=True)
    
    dest_g = GLOBAL_SKILLS / skill_id
    dest_w = WS_SKILLS / skill_id
    
    if dest_g.exists():
        shutil.rmtree(dest_g)
    shutil.copytree(src, dest_g)
    
    if dest_w.exists():
        shutil.rmtree(dest_w)
    shutil.copytree(src, dest_w)
    
    print(f"Successfully installed skill '{skill_id}' to:")
    print(f"  - Global: {dest_g}")
    print(f"  - Workspace: {dest_w}")

def show_info(skill_id):
    skill_file = SKILLS_DIR / skill_id / "SKILL.md"
    if not skill_file.exists():
        print(f"Skill file not found: {skill_file}")
        return
    print(skill_file.read_text(encoding="utf-8", errors="replace"))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(0)
    
    action = sys.argv[1].lower()
    arg = sys.argv[2]
    
    if action == "search":
        search_skills(arg)
    elif action == "install":
        install_skill(arg)
    elif action == "info":
        show_info(arg)
    else:
        print(f"Unknown action: {action}")
