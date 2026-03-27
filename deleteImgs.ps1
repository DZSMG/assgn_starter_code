# Remove the 'images' directory if it exists
if (Test-Path "images") {
    Remove-Item -Recurse -Force "images"
}

# Create the 'images' directory
New-Item -ItemType Directory -Path "images" | Out-Null

# Remove the 'result.mp4' file if it exists
if (Test-Path "result.mp4") {
    Remove-Item -Force "result.mp4"
}