import sys
import pyfiglet
from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich.align import Align

# Bash: sudo pip3 install rich pyfiglet


print("")

console = Console()

# Generate ASCII Art for the name
ascii_art = pyfiglet.figlet_format("WES Pipeline", font="slant")

# Create the content
content = Text(ascii_art, style="bold blue")
content.append("WES Analysis Pipeline ", style="bold white")
content.append("(v1.0.0 | Build 2025)\n", style="dim white")

# Metadata section
meta_text = Text()
meta_text.append("Made with ❤️  by Muhammad Ashraf\n", style="bold italic white")
meta_text.append("GitHub Repo: ", style="bold white")
meta_text.append("github.com/muhammadash717/WES-Pipeline\n", style="underline blue")

# Combine
full_content = Text.assemble(content, meta_text)

# Create a panel with a DNA-themed border
panel = Panel(
    Align.center(full_content),
    border_style="bold blue",    
    subtitle=f"[bold white]Sample : {sys.argv[1]}" if len(sys.argv) > 1 else "",
    expand = False, padding = (0, 5)
)

console.print(panel)

print("")