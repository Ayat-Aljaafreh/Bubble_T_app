# main_driver.py

# استدعاء ملفات bubble و dew
import bubble_T
import dew_T     

def main():
    print("Choose calculation type:")
    print("1 - Bubble Point")
    print("2 - Dew Point")

    while True:
        choice = input("Enter 1 or 2: ").strip()
        if choice == "1":
            bubble_T.main()   # يشغّل ملف bubble
            break
        elif choice == "2":
            dew_T.main()      # يشغّل ملف dew
            break
        else:
            print("Invalid choice. Enter 1 for Bubble or 2 for Dew.")

if __name__ == "__main__":
    main()
