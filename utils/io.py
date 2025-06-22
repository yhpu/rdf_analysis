def get_user_selection(prompt, options):
    print(f"{prompt} from: {options}")
    while True:
        choice = input(">> ").strip()
        if choice in options:
            return choice
        print("Invalid choice. Try again.")


def yes_or_no(prompt):
    while True:
        ans = input(f"{prompt} (yes/no): ").strip().lower()
        if ans in ("yes", "y"):
            return True
        elif ans in ("no", "n"):
            return False
        print("Please enter yes or no.")