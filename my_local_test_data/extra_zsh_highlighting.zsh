# Extra zsh-syntax-highlighting styles and helpers
# -------------------------------------------------
#  This file lives inside the project workspace so that **no user dot-files are
#  modified automatically**.  To activate it, add a single line to your own
#  ~/.zshrc (or ~/.z-monokai) **after** the line that sources
#  zsh-syntax-highlighting:
#
#     source /ABS/PATH/TO/barskilab-workflows/my_local_test_data/extra_zsh_highlighting.zsh
#
#  Then restart the terminal or run:   source ~/.zshrc
# -------------------------------------------------

# Ensure the helper highlighters are enabled.  Keep 'main' first.
if [[ -z "${ZSH_HIGHLIGHT_HIGHLIGHTERS+x}" ]]; then
  ZSH_HIGHLIGHT_HIGHLIGHTERS=(main brackets cursor root pattern)
else
  # Unify the array and append missing highlighters.
  typeset -U ZSH_HIGHLIGHT_HIGHLIGHTERS
  ZSH_HIGHLIGHT_HIGHLIGHTERS+=(brackets cursor root pattern)
fi

# Initialize pattern highlighting (force enable)
ZSH_HIGHLIGHT_PATTERNS=()
# Highlight TODO / FIXME inside comments
ZSH_HIGHLIGHT_PATTERNS+=('TODO|FIXME|BUG' 'fg=red,bold')

# -------------------------------------------------
# Professional Muted Palette - Optimized for eye comfort and readability
# -------------------------------------------------

# ===========================================
# COMMANDS (muted teal) - "things you execute"
# High frequency → comfortable, readable color
# ===========================================
ZSH_HIGHLIGHT_STYLES[command]='fg=74'                     # External commands (teal)
ZSH_HIGHLIGHT_STYLES[builtin]='fg=74,bold'                # Shell builtins (cd, echo)
ZSH_HIGHLIGHT_STYLES[function]='fg=80'                    # User functions (lighter teal)
ZSH_HIGHLIGHT_STYLES[alias]='fg=74'                       # Aliases (ls -> lsd)
ZSH_HIGHLIGHT_STYLES[suffix-alias]='fg=74,underline'      # Suffix aliases
ZSH_HIGHLIGHT_STYLES[global-alias]='fg=74,bold'           # Global aliases

# ===========================================
# OPERATORS (muted purple) - "how commands connect/flow"
# High frequency → subtle but visible
# ===========================================
ZSH_HIGHLIGHT_STYLES[commandseparator]='fg=141'           # && || | ; (muted purple)
ZSH_HIGHLIGHT_STYLES[redirection]='fg=141,bold'           # > >> < <<<
ZSH_HIGHLIGHT_STYLES[process-substitution]='fg=141'       # <(command)
ZSH_HIGHLIGHT_STYLES[process-substitution-delimiter]='fg=141,bold'

# ===========================================
# ARGUMENTS (warm orange) - "how commands behave"
# High frequency → replaces harsh yellow with readable orange
# ===========================================
ZSH_HIGHLIGHT_STYLES[single-hyphen-option]='fg=208'       # -a -v (warm orange)
ZSH_HIGHLIGHT_STYLES[double-hyphen-option]='fg=208,bold'  # --help --version
ZSH_HIGHLIGHT_STYLES[assign]='fg=208'                     # VAR=value

# ===========================================
# DATA (green family) - "what commands work with"
# Distinguished by functionality: literal vs dynamic
# ===========================================
ZSH_HIGHLIGHT_STYLES[single-quoted-argument]='fg=108'     # 'string' (dim green - literal)
ZSH_HIGHLIGHT_STYLES[double-quoted-argument]='fg=114,bold' # "string" (bright green - dynamic)
ZSH_HIGHLIGHT_STYLES[dollar-quoted-argument]='fg=114'     # $'string' (medium green)

# ===========================================
# VARIABLES (bright blue) - "dynamic data"
# Medium frequency → needs to stand out more
# ===========================================
ZSH_HIGHLIGHT_STYLES[dollar-double-quoted-argument]='fg=39' # $VAR in quotes (bright blue)
ZSH_HIGHLIGHT_STYLES[back-quoted-argument]='fg=39,bold'     # `command`
ZSH_HIGHLIGHT_STYLES[command-substitution]='fg=39'         # $(command)
ZSH_HIGHLIGHT_STYLES[command-substitution-delimiter]='fg=39,bold'
ZSH_HIGHLIGHT_STYLES[arithmetic-expansion]='fg=39'         # $((2+2))
ZSH_HIGHLIGHT_STYLES[history-expansion]='fg=39,bold'       # !! !$

# Bare variables (the key fix for $HOME)
ZSH_HIGHLIGHT_STYLES[default]='fg=39'                     # $VAR standalone variables

# ===========================================
# PATHS (light gray) - "file system references"
# Low frequency → subtle but distinct from default text
# ===========================================
ZSH_HIGHLIGHT_STYLES[path]='fg=250,underline'             # Light gray
ZSH_HIGHLIGHT_STYLES[path_pathseparator]='fg=250,bold'
ZSH_HIGHLIGHT_STYLES[path_prefix]='fg=250'
ZSH_HIGHLIGHT_STYLES[path_prefix_pathseparator]='fg=250,bold'

# ===========================================
# BRACKETS (muted rainbow) - "structure indicators"
# Low frequency → can use brighter colors for distinction
# ===========================================
ZSH_HIGHLIGHT_STYLES[bracket-level-1]='fg=74,bold'        # Teal
ZSH_HIGHLIGHT_STYLES[bracket-level-2]='fg=114,bold'       # Green
ZSH_HIGHLIGHT_STYLES[bracket-level-3]='fg=141,bold'       # Purple
ZSH_HIGHLIGHT_STYLES[bracket-level-4]='fg=208,bold'       # Orange
ZSH_HIGHLIGHT_STYLES[cursor-matchingbracket]='fg=255,bold,bg=39'

# ===========================================
# DANGEROUS COMMANDS - "elevated privileges & destructive operations"
# ===========================================
ZSH_HIGHLIGHT_STYLES[precommand]='fg=208,bg=52'            # sudo (orange on dark red bg)

# ===========================================
# NESTED QUOTES - "improved string parsing"
# ===========================================
ZSH_HIGHLIGHT_STYLES[single-quoted-argument-unclosed]='fg=114,bg=52'    # Unclosed quotes
ZSH_HIGHLIGHT_STYLES[double-quoted-argument-unclosed]='fg=114,bg=52'    # Unclosed quotes
ZSH_HIGHLIGHT_STYLES[dollar-quoted-argument-unclosed]='fg=114,bg=52'    # Unclosed quotes

# ===========================================
# ERRORS (red family) - "problems"
# ===========================================
ZSH_HIGHLIGHT_STYLES[bracket-error]='fg=red,bold'
ZSH_HIGHLIGHT_STYLES[unknown-token]='fg=red,bold'

# ===========================================
# PATTERN-BASED ENHANCEMENTS - Advanced semantic highlighting
# ===========================================
# Dangerous commands - visual warning (orange on dark red)
ZSH_HIGHLIGHT_PATTERNS+=('\brm\b' 'fg=208,bg=52')
ZSH_HIGHLIGHT_PATTERNS+=('\bmv\b' 'fg=208,bg=52') 
ZSH_HIGHLIGHT_PATTERNS+=('\bchmod\b' 'fg=208,bg=52')
ZSH_HIGHLIGHT_PATTERNS+=('\bchown\b' 'fg=208,bg=52')
ZSH_HIGHLIGHT_PATTERNS+=('\bdd\b' 'fg=208,bg=52')

# Network commands - bright cyan
ZSH_HIGHLIGHT_PATTERNS+=('\bssh\b' 'fg=81')
ZSH_HIGHLIGHT_PATTERNS+=('\bcurl\b' 'fg=81')
ZSH_HIGHLIGHT_PATTERNS+=('\bwget\b' 'fg=81')
ZSH_HIGHLIGHT_PATTERNS+=('\brsync\b' 'fg=81')
ZSH_HIGHLIGHT_PATTERNS+=('\bscp\b' 'fg=81')

# Escape sequences in strings
ZSH_HIGHLIGHT_PATTERNS+=('\\[ntrfbav]' 'fg=39,bold')

# Security context awareness
ZSH_HIGHLIGHT_PATTERNS+=('\bsudo\b' 'fg=208,bg=52,bold')

# Dangerous combinations
ZSH_HIGHLIGHT_PATTERNS+=('\brm\s+.*-rf\b' 'fg=15,bg=196,bold')  # rm -rf = white on bright red

# ===========================================
# NESTED QUOTES HIGHLIGHTING - Advanced string parsing
# ===========================================
# Single quotes inside double quotes: "text 'nested' text"
ZSH_HIGHLIGHT_PATTERNS+=('"[^"]*'\''[^'\'']*'\''[^"]*"' 'fg=114,bold,bg=22')  # Green with dark green bg

# Double quotes inside single quotes: 'text "nested" text'  
ZSH_HIGHLIGHT_PATTERNS+=('''[^'\'']*"[^"]*"[^'\'']*''' 'fg=108,bold,bg=58')   # Dim green with brown bg

# Multiple nested quotes: "outer 'inner "deep" inner' outer"
ZSH_HIGHLIGHT_PATTERNS+=('"[^"]*'\''[^'\'']*"[^"]*"[^'\'']*'\''[^"]*"' 'fg=15,bg=21,bold')  # White on dark blue

# Escaped quotes within strings
ZSH_HIGHLIGHT_PATTERNS+=('\\"' 'fg=226,bold')  # Escaped double quote - bright yellow
ZSH_HIGHLIGHT_PATTERNS+=('''\\'\''' 'fg=226,bold')  # Escaped single quote - bright yellow

# ===========================================
# ADVANCED PRO-LEVEL FEATURES
# ===========================================
# Context-aware highlighting based on current directory
if [[ "$PWD" =~ "^/(etc|usr|var|sys)" ]]; then
  # In system directories - highlight more aggressively
  ZSH_HIGHLIGHT_PATTERNS+=('\b(cp|mv|rm|chmod|chown)\b' 'fg=15,bg=196,bold')
fi

# Multi-line command intelligence
ZSH_HIGHLIGHT_PATTERNS+=('\\\s*$' 'fg=39,bold')  # Line continuation backslash

# Error prediction patterns
ZSH_HIGHLIGHT_PATTERNS+=('\$\{[A-Z_]+\}' 'fg=196,bold')  # Undefined variables pattern
ZSH_HIGHLIGHT_PATTERNS+=('\brm\s+.*\$' 'fg=15,bg=196,bold')  # rm with variable at end

# File extension awareness in arguments
ZSH_HIGHLIGHT_PATTERNS+=('\b\S+\.(sh|py|js|rb|go|rs)\b' 'fg=81,underline')  # Script files
ZSH_HIGHLIGHT_PATTERNS+=('\b\S+\.(txt|md|doc|pdf)\b' 'fg=250,underline')   # Document files
ZSH_HIGHLIGHT_PATTERNS+=('\b\S+\.(jpg|png|gif|svg)\b' 'fg=206,underline')  # Image files

# Advanced operators
ZSH_HIGHLIGHT_PATTERNS+=('\b(xargs|parallel|find)\b' 'fg=81,bold')  # Pipeline power tools

# Git-specific highlighting
if [[ -d ".git" ]] || git rev-parse --git-dir > /dev/null 2>&1; then
  ZSH_HIGHLIGHT_PATTERNS+=('\bgit\s+(reset|rebase|push\s+--force)\b' 'fg=208,bg=52,bold')
fi